#!/usr/bin/env python3
"""Sideratlas v1.2: restrained dynamic accretion and a neutral geodesic plunge.

Null paths, emitter motion, frequency shifts, delays and capture are calculated.
Brightness structure, plunge heating and its optical depth are explicit source
prescriptions. No screen-space rings, bloom or image-order boosts are added.
The v1.1 renderer remains available unchanged in render_kn_orbit.py.
"""
from __future__ import annotations
import argparse
import gc
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
from PIL import Image, ImageDraw

from render_kn_orbit import (KerrNewmanRenderer, view_angle_at_phase, ROOT, FPS,
    FOV, SPIN, CHARGE, OUTER_RADIUS, INCLINATION, LOOP_TG, EXPOSURE,
    DEFAULT_SUPERSAMPLE, LOOP_DURATION, REPEATS, DURATION, VIDEO_NAMES, LOOP_VIDEO_NAMES)
from kn_accretion_trace import KerrNewmanAccretionTracer
from kn_flow import NeutralAccretionFlow
from kn_dynamic_emission import DynamicDiskEmission
from instagram_export import encoder

OUTPUT_DIR = ROOT/'output-sideratlas-kn-v120'
PLUNGE_OPTICAL_DEPTH = .15
PHYSICAL_SOURCES = ('render_kn_accretion.py','kn_accretion_trace.py',
    'kn_accretion_trace.cpp','kn_flow.py','kn_dynamic_emission.py',
    'kn_orbits.py','kn_orbit_emission.py','kn_orbit_trace.py','nasa_emission.py')


def source_hash():
    digest = hashlib.sha256()
    for name in (*PHYSICAL_SOURCES,'render_kn_orbit.py','instagram_export.py',
                 'render.py','emission.py','raytrace.py'):
        digest.update((ROOT/name).read_bytes())
    for path in sorted((ROOT/'fonts').glob('*.otf')):
        digest.update(path.read_bytes())
    return digest.hexdigest()


def compact_hits(maps, exclude=None):
    """Discard unused event slots before flow/emission preparation."""
    mask = maps['mask'].copy()
    if exclude is not None:
        mask &= ~np.asarray(exclude)[...,None]
    index = np.flatnonzero(mask.ravel())
    slots = mask.shape[-1]
    pixel = index//slots
    result = {'pixel':pixel.astype(np.int32),'slot':(index%slots).astype(np.uint8)}
    for key in ('r','phi','delay_tg','kr_over_E','hit_layer','hit_kind'):
        result[key] = maps[key].ravel()[index]
    for key in ('lambda','eta'):
        result[key] = maps[key].ravel()[pixel]
    return result


def trace_resolved(tracer, alpha, beta, angle):
    """Retry only rare unresolved rays; never silently cap an image order."""
    maps = tracer.trace_samples(alpha,beta,angle,allow_unresolved=True)
    bad = maps['status'] >= 3
    events = compact_hits(maps,bad)
    stats = {'rays':int(bad.size),'retried_rays':int(bad.sum()),
        'maximum_constraint_residual':float(maps['constraint_residual'].max()),
        'captured_rays':int(np.count_nonzero(maps['status']==2))}
    if bad.any():
        indices = np.flatnonzero(bad)
        retry = KerrNewmanAccretionTracer(8,8,spin=tracer.spin,charge_q=tracer.charge_q,
            rin=tracer.rin,rout=tracer.rout,angle_step=min(.006,tracer.angle_step/2),
            max_angle=40*math.pi,max_hits=16,horizon_offset_rg=tracer.horizon_offset_rg,
            threads=tracer.threads)
        retraced = retry.trace_samples(np.asarray(alpha).ravel()[indices],
                                      np.asarray(beta).ravel()[indices],angle)
        extra = compact_hits(retraced)
        extra['pixel'] = indices[extra['pixel']].astype(np.int32)
        for key in events:
            events[key] = np.concatenate((events[key],extra[key]))
        stats['maximum_constraint_residual'] = max(stats['maximum_constraint_residual'],
            float(retraced['constraint_residual'].max()))
        stats['captured_rays'] += int(np.count_nonzero(retraced['status']==2))
    stats['maximum_image_order'] = int(events['hit_layer'].max()) if len(events['r']) else -1
    stats['image_order_hits'] = {str(int(order)):int(count) for order,count in zip(
        *np.unique(events['hit_layer'],return_counts=True))}
    stats['plunging_events'] = int(np.count_nonzero(events['hit_kind']==2))
    del maps
    return events,stats


def transfer_weights(events, flow, pixel_count, optical_depth=PLUNGE_OPTICAL_DEPTH):
    """Ordered homogeneous thin-slab transfer, followed by g^4 intensity.

    The supplied plunge intensity is its face-on emergent intensity. Its
    source function is I_face/(1-exp(-tau)); tau is an illustrative normal
    optical depth. The opaque outer surface has isotropic Page--Thorne I.
    """
    g = flow.redshift(events['r'],events['lambda'],events['kr_over_E'])
    mu = g*np.sqrt(np.maximum(events['eta'],0))/events['r']
    if np.any(mu > 1+3e-5) or not np.isfinite(g).all():
        raise AssertionError('invalid emitter-frame photon direction or shift')
    mu = np.clip(mu,0,1)
    inner = events['hit_kind']==2
    depth = np.zeros_like(g)
    # Once tau/mu exceeds 80, residual transmission is <2e-35; this is a
    # numerical exponential bound, not a grazing-angle brightness floor.
    depth[inner] = np.minimum(80,optical_depth/np.maximum(mu[inner],1e-300))
    transmission = np.where(inner,np.exp(-depth),0.)
    emitted = np.where(inner,-np.expm1(-depth)/(-math.expm1(-optical_depth)),1.)
    foreground = np.ones(pixel_count,dtype=np.float64)
    weights = np.zeros_like(g)
    for slot in np.unique(events['slot']):
        select = events['slot']==slot
        pix = events['pixel'][select]
        weights[select] = foreground[pix]*emitted[select]*g[select]**4
        foreground[pix] *= transmission[select]
    if not np.isfinite(weights).all() or np.any(weights<0):
        raise AssertionError('nonfinite transfer weights')
    return weights, g


def critical_sampling_mask(width,height,fov,angle,spin=SPIN,charge=CHARGE):
    """Use the spherical-null critical curve solely to allocate extra rays.

    R=R'=0 gives P=4rDelta/Delta', lambda=(r^2+a^2-P)/a,
    eta=16r^2Delta/Delta'^2-(lambda-a)^2. No curve is painted into RGB.
    """
    from kn_orbits import horizon_radius
    image = Image.new('1',(width,height)); draw=ImageDraw.Draw(image)
    theta = math.radians(angle); si,ci = math.sin(theta),math.cos(theta)
    rp = horizon_radius(spin,charge)
    if abs(spin)<1e-12:
        photon_radius=(3+math.sqrt(9-8*charge*charge))/2
        delta=photon_radius**2-2*photon_radius+charge*charge
        radius=photon_radius**2/math.sqrt(delta)
        phi=np.linspace(0,2*math.pi,2049)
        draw.line(list(zip(width/2+radius*np.cos(phi)*width/fov,
                           height/2-radius*np.sin(phi)*width/fov)),fill=1,width=7)
        return np.asarray(image,dtype=bool).copy()
    rr = np.linspace(rp+1e-6,8.,32768)
    delta = rr*rr-2*rr+spin*spin+charge*charge
    P = 2*rr*delta/(rr-1)
    lam = (rr*rr+spin*spin-P)/spin
    K = 4*rr*rr*delta/(rr-1)**2
    eta = K-(lam-spin)**2
    if abs(si)<.03:
        # A polar observer has lambda=0. The small-angle circular seed only
        # locates a conservative sampling band; the traced rays remain exact.
        roots = np.flatnonzero(lam[:-1]*lam[1:]<=0)
        if not roots.size:
            raise RuntimeError('failed to locate polar spherical photon orbit')
        j=int(roots[0]); ll,rad=lam[j:j+2],rr[j:j+2]
        if ll[0]>ll[1]:ll,rad=ll[::-1],rad[::-1]
        r=float(np.interp(0,ll,rad))
        dd=r*r-2*r+spin*spin+charge*charge
        radius=math.sqrt(4*r*r*dd/(r-1)**2)
        phi=np.linspace(0,2*math.pi,2049)
        xs,ys=radius*np.cos(phi),radius*np.sin(phi)
        points=list(zip(width/2+xs*width/fov,height/2-ys*width/fov))
        draw.line(points,fill=1,width=9)
    else:
        beta2=eta+spin*spin*ci*ci-lam*lam*ci*ci/(si*si)
        valid=(beta2>=0)&(eta>=-spin*spin)&np.isfinite(beta2)
        alpha=-lam[valid]/si; beta=np.sqrt(beta2[valid])
        for sign in (1,-1):
            points=list(zip(width/2+alpha*width/fov,height/2-sign*beta*width/fov))
            if len(points)>1:draw.line(points,fill=1,width=7)
    return np.asarray(image,dtype=bool).copy()


class PreparedRayEmission:
    def __init__(self,events,source,flow,pixel_count):
        self.events,self.source,self.pixel_count=events,source,pixel_count
        self.weights,self.g=transfer_weights(events,flow,pixel_count)
        self.prepared=source.prepare(events['r'],events['phi'],events['delay_tg'])
        self.opposite_prepared=None

    def evaluate(self,time_tg,opposite=False,structured=True):
        if opposite and self.opposite_prepared is None:
            self.opposite_prepared=self.source.prepare(self.events['r'],
                np.remainder(self.events['phi'].astype(float)+math.pi,2*math.pi),
                self.events['delay_tg'])
        prepared=self.opposite_prepared if opposite else self.prepared
        intensities=prepared.evaluate(time_tg,structured=structured)
        return np.bincount(self.events['pixel'],weights=intensities*self.weights,
                           minlength=self.pixel_count)


class KerrNewmanAccretionRenderer(KerrNewmanRenderer):
    def __init__(self,width=1080,supersample=DEFAULT_SUPERSAMPLE,
                 inclination=INCLINATION,fov=FOV,spin=SPIN,charge=CHARGE,
                 angle_step=.012,refinement=2,temporal_samples=3):
        if width<270 or width%2 or supersample<1 or refinement not in (1,2,3,4):
            raise ValueError('even width>=270, sampling>=1 and refinement1..4 required')
        if temporal_samples not in (1,3):raise ValueError('temporal samples must be1or3')
        self.width,self.height=width,round(width*16/9)//2*2
        self.scale,self.fov,self.spin,self.charge=width/1080,float(fov),float(spin),float(charge)
        self.inclination,self.refinement,self.temporal_samples=inclination,refinement,temporal_samples
        self._configure_layout()
        self.flow=NeutralAccretionFlow(spin_a_star=spin,charge_q=charge)
        self.source=DynamicDiskEmission(flow=self.flow,spin_a_star=spin,charge_q=charge)
        self.tracer=KerrNewmanAccretionTracer(round(width*supersample),round(self.hero_h*supersample),
            fov=fov,spin=spin,charge_q=charge,angle_step=angle_step,max_hits=4)
        self.refined_tracer=KerrNewmanAccretionTracer(8,8,spin=spin,charge_q=charge,
            angle_step=min(.006,angle_step),max_hits=6)
        self.geometry={k:v for k,v in self.tracer.metadata.items() if not k.startswith('last_')}
        self.current_inclination=None
        self.overlays={mode:self._overlay(mode) for mode in VIDEO_NAMES}
        self.set_view(0)

    def set_view(self,phase):
        angle=view_angle_at_phase(phase,self.inclination)
        if self.current_inclination is not None and abs(angle-self.current_inclination)<1e-12:return
        # Release old hit-only source preparations before allocating another
        # dense event map; large unused max_hits slots never enter emission.
        self.coarse=self.fine=None
        gc.collect()
        events,stats=trace_resolved(self.tracer,self.tracer.x,self.tracer.y,angle)
        h,w=self.tracer.height,self.tracer.width
        self.coarse=PreparedRayEmission(events,self.source,self.flow,h*w)
        refine=critical_sampling_mask(w,h,self.fov,angle,self.spin,self.charge)
        important=(events['hit_layer']>=1)|(events['hit_kind']==2)
        refine.ravel()[events['pixel'][important]]=True
        self.refined_pixels=np.flatnonzero(refine) if self.refinement>1 else np.array([],int)
        if len(self.refined_pixels):
            n=self.refinement
            offsets=(np.arange(n)+.5)/n-.5
            ox,oy=np.meshgrid(offsets,offsets)
            spacing=self.fov/w
            alpha=self.tracer.x.ravel()[self.refined_pixels,None]+ox.ravel()*spacing
            beta=self.tracer.y.ravel()[self.refined_pixels,None]+oy.ravel()*spacing
            finer,fine_stats=trace_resolved(self.refined_tracer,alpha,beta,angle)
            self.fine=PreparedRayEmission(finer,self.source,self.flow,alpha.size)
            stats['refined']=fine_stats
        stats['refined_parent_pixels']=len(self.refined_pixels)
        stats['refinement_grid']=self.refinement
        stats['base_observed_g_range']=[float(self.coarse.g.min()),float(self.coarse.g.max())]
        self.view_diagnostics=stats
        self.current_inclination=angle

    def _radiance_at_time(self,time_tg,opposite=False,structured=True):
        values=self.coarse.evaluate(time_tg,opposite,structured)
        if self.fine is not None:
            finer=self.fine.evaluate(time_tg,opposite,structured)
            values[self.refined_pixels]=finer.reshape(-1,self.refinement**2).mean(axis=1)
        image=values.reshape(self.tracer.height,self.tracer.width)
        return np.fliplr(image) if opposite else image

    def _exposed_radiance(self,phase,opposite=False,structured=True):
        t=(phase%1)*LOOP_TG
        if self.temporal_samples==1:return self._radiance_at_time(t,opposite,structured).astype('float32')
        # Three-point Gaussian quadrature of a 180-degree source shutter.
        # Central-view geodesics are held fixed: no claim of camera-motion blur.
        half_width=LOOP_TG/(FPS*LOOP_DURATION)*.25
        offset=math.sqrt(3/5)*half_width
        result=np.zeros((self.tracer.height,self.tracer.width),float)
        for dt,weight in ((-offset,5/18),(0,4/9),(offset,5/18)):
            result+=weight*self._radiance_at_time(t+dt,opposite,structured)
        if not np.isfinite(result).all() or np.any(result<0):raise AssertionError('invalid exposed radiance')
        return result.astype('float32')

    def radiance(self,phase,structured=True):
        self.set_view(phase)
        return self._exposed_radiance(phase,False,structured)

    def antipodal_radiance(self,phase):
        self.set_view(phase)
        return self._exposed_radiance(phase+.5,True)

    def metadata(self):
        self.set_view(0)
        return {'version':'1.2.0','width':self.width,'height':self.height,'fps':FPS,
            'frames':FPS*DURATION,'duration_seconds':DURATION,'loop_duration_seconds':LOOP_DURATION,
            'repetitions':REPEATS,'unique_frames':FPS*LOOP_DURATION,
            'view_angle_degrees':[self.inclination-360,self.inclination],
            'camera_orbit':{'law':'theta=83-360*phase; continuous oriented meridional scan',
                'angular_span_deg':-360,'observer':'sequence of static observers at infinity; no camera-velocity aberration'},
            'spin_configuration':{'a_star':self.spin,'charge_q':self.charge,
                'alignment':'aligned, prograde equatorial neutral test matter',
                'choice':'hypothetical large charge, not a fit to a typical astrophysical black hole'},
            'geometry':self.geometry,'flow':self.flow.metadata(),'emission':self.source.metadata(),
            'radiative_transfer':{'bolometric':'g^4 I_em; g^-1=ut-lambda*uphi-(k_r/E)*ur',
                'time':'t_emit=t_observer-delay independently at every intersection',
                'outer_disk':'first opaque isotropic Page-Thorne surface',
                'plunge_normal_optical_depth':PLUNGE_OPTICAL_DEPTH,
                'plunge_slab':'homogeneous thin slab; I_em(mu)=I_face*(1-exp(-tau/mu))/(1-exp(-tau)); foreground transmission exp(-tau/mu)',
                'opacity_is_prescribed':True,'horizon':'no source or surface on or inside r_plus; captured rays can see foreground gas'},
            'sampling':{'refinement_grid':self.refinement,
                'refinement_selection':'spherical-null critical curve sampling band plus every higher-order or plunging-gas parent pixel',
                'refinement':'independently traced subpixel null geodesics; equal-area average in radiance; no painted rings',
                'overflow':'sparse retrace with16eventslots and half RK4 step; unresolved rays raise',
                'temporal_samples':self.temporal_samples,
                'shutter':'180-degree source-only exposure,3pointGauss quadrature; fixed central-view photon geometry'},
            'layout':{'brand':'@sideratlas','full_frame':True,'hero_x':0,'hero_y':0,
                'hero_width':self.width,'hero_height':self.height,'header_top':270,
                'display_roll_degrees':self.display_roll,'viewport_fov_horizontal_rg':self.viewport_fov,
                'lower_text_in_caption':True,'caption_file':'caption-kerr-newman.text'},
            'display':{'exposure':EXPOSURE,'curve':'x/(1+x), x=radiance/exposure',
                'color':'illustrative orange bolometric brightness, not a visible spectrum',
                'averaging':'native radiance plus selective independent subpixel rays; rotated viewport filtered before display curve'},
            'disk_faces':'identical equatorial source and transfer on both faces',
            'geodesic_reuse':'120independent central views; exact antipodal geometry with new source evaluation at phi+pi and own retardedtime',
            'time':'800t_g per8-second constructed source/camera loop, four identical encoded repeats',
            'outputs':VIDEO_NAMES,'loop_outputs':LOOP_VIDEO_NAMES,
            'encoding':__import__('instagram_export').metadata(),
            'source_sha256':source_hash()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width',type=int,default=1080)
    parser.add_argument('--supersample',type=float,default=DEFAULT_SUPERSAMPLE)
    parser.add_argument('--refinement',type=int,default=2)
    parser.add_argument('--temporal-samples',type=int,default=3)
    parser.add_argument('--preview',action='store_true')
    parser.add_argument('--output',type=Path,default=OUTPUT_DIR)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    print('Tracing dynamic accretion, neutral plunge and multiple Kerr-Newman images',flush=True)
    renderer=KerrNewmanAccretionRenderer(args.width,args.supersample,refinement=args.refinement,
                                       temporal_samples=args.temporal_samples)
    first=renderer.frames(0)
    for mode,frame in first.items():frame.save(args.output/f'sideratlas-kerr-newman-{mode}-cover.png')
    renderer.frames(0,structured=False)['color'].save(args.output/'steady-reference.png')
    sheet=Image.new('RGB',(1080,1280))
    for col,phase in enumerate((0,.25,.5)):
        frames=first if phase==0 else renderer.frames(phase)
        for row,mode in enumerate(VIDEO_NAMES):
            sheet.paste(frames[mode].resize((360,640),Image.Resampling.LANCZOS),(col*360,row*640))
    sheet.save(args.output/'contact-sheet.jpg',quality=95)
    meta=renderer.metadata();meta['supersample']=args.supersample
    signature=hashlib.sha256(json.dumps(meta,sort_keys=True).encode()).hexdigest()
    meta['radiance_signature']=signature
    (args.output/'render-metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    if args.preview:
        (args.output/'view-diagnostics.json').write_text(json.dumps(renderer.view_diagnostics,indent=2)+'\n')
        print(f'Preview finished in{time.perf_counter()-started:.1f}s',flush=True);return
    cache=args.output/'.render-cache'/signature[:20];cache.mkdir(parents=True,exist_ok=True)
    total=FPS*LOOP_DURATION
    processes={mode:encoder(args.output/(Path(name).stem+'.partial.mp4'),renderer.width,renderer.height,FPS)
               for mode,name in LOOP_VIDEO_NAMES.items()}
    diagnostic_path=cache/'view-diagnostics.json'
    diagnostics=json.loads(diagnostic_path.read_text()) if diagnostic_path.exists() else []
    try:
        for n in range(total):
            path=cache/f'{n:03d}.npy'
            if path.exists():
                values=np.load(path,allow_pickle=False)
                if values.shape!=(renderer.tracer.height,renderer.tracer.width) or not np.isfinite(values).all() or np.any(values<0):raise ValueError('invalid resume frame')
            else:
                values=renderer.radiance(n/total)
                tmp=path.with_suffix('.partial.npy');np.save(tmp,values,allow_pickle=False);tmp.replace(path)
                if n<total//2:
                    opposite=cache/f'{n+total//2:03d}.npy'
                    if not opposite.exists():
                        other=renderer.antipodal_radiance(n/total)
                        tmp=opposite.with_suffix('.partial.npy');np.save(tmp,other,allow_pickle=False);tmp.replace(opposite)
                    diagnostics=[entry for entry in diagnostics if entry['frame']!=n]
                    diagnostics.append({'frame':n,**renderer.view_diagnostics})
                    diagnostics.sort(key=lambda entry:entry['frame'])
                    diagnostic_path.write_text(json.dumps(diagnostics,indent=2)+'\n')
                    (args.output/'view-diagnostics.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
            for mode,process in processes.items():process.stdin.write(renderer.compose(values,mode).tobytes())
            if n%15==0:print(f'Frame{n+1}/{total};elapsed{time.perf_counter()-started:.0f}s',flush=True)
    finally:
        for process in processes.values():process.stdin.close()
        codes=[process.wait() for process in processes.values()]
    if any(codes):raise RuntimeError(f'incomplete encodes:{codes}')
    for name in LOOP_VIDEO_NAMES.values():(args.output/(Path(name).stem+'.partial.mp4')).replace(args.output/name)
    for mode,name in VIDEO_NAMES.items():
        subprocess.run([shutil.which('ffmpeg') or 'ffmpeg','-v','error','-y','-stream_loop',str(REPEATS-1),
            '-i',str(args.output/LOOP_VIDEO_NAMES[mode]),'-map','0:v:0','-an','-c:v','copy',
            '-movflags','+faststart','-use_editlist','0',str(args.output/name)],check=True)
    manifest={'radiance_signature':signature,'source_render_sha256':meta['source_sha256'],
        'frames_sha256':{f'{n:03d}.npy':hashlib.sha256((cache/f'{n:03d}.npy').read_bytes()).hexdigest() for n in range(total)}}
    (cache/'radiance-cache-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (args.output/'view-diagnostics.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    meta['elapsed_seconds']=time.perf_counter()-started
    (args.output/'render-metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(f'Finished in{meta["elapsed_seconds"]:.1f}s',flush=True)


if __name__=='__main__':main()
