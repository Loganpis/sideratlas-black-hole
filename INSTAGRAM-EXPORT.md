# Instagram delivery for @sideratlas

The delivery file is a native **1080 × 1920, 30 fps, 32-second SDR MP4** (four copies of one eight-second cycle). It is encoded directly from the animation's RGB frames, without an intermediate compressed video. The disk is traced at approximately the final angular pixel density; the rotated viewport is filtered in floating-point radiance before display encoding. The export does not upscale the finished video, sharpen the rings, or add grain.

Version 1.1.0 uses the entire portrait frame for a close diagonal view. The default source detector is 2700 × 2250, covering an unrolled field of view of 72 × 60 r_g. Floating-point radiance is bilinearly sampled into a 2× portrait viewport with a 28 r_g horizontal field of view and a fixed −60.642° screen rotation, then box-averaged to delivery size before tone mapping. The close framing deliberately crops some outer-disk emission in certain views, especially near the poles. This is a presentation change; the physical model is unchanged.

Only **@sideratlas**, **NO.001**, and **THE BLACK HOLE** remain on screen, over the simulation, with a subtle glyph shadow. All former lower text, model disclosures, and source credits are in [`caption-kerr-newman.text`](caption-kerr-newman.text), with a compatible `.txt` copy. Use that description when posting the video so the scientific context accompanies the clean visual.

## Export choices

`instagram_export.py` uses H.264 High, 4:2:0 chroma, the slow preset, CRF 14, a 20 Mbps maximum rate and a 40 Mbit buffer. It writes a constant frame rate, a closed two-second GOP with no B-frames, a front-loaded MP4 movie header, and no edit lists. The video is silent. CRF 14 and these rate-control choices are conservative production decisions for fine rings and typography; they are **not** an Instagram quality requirement or a promise that the platform will retain these settings. The average bitrate can be much lower than its cap, particularly with a black background and static text.

Input frames use full-range sRGB values. FFmpeg's zscale filter performs a real sRGB-to-Rec.709 transfer conversion and RGB-to-limited-range YCbCr conversion. This is a display-referred conversion: zimg uses the sRGB display EOTF and the inverse BT.1886 gamma-2.4 display curve for Rec.709, preserving the intended display light of the rendered graphics. It is not a camera/scene OETF conversion. Output primaries, transfer and matrix are all tagged BT.709, with limited-range video and left-sited chroma. Expanding the RGB input to 16-bit planar storage avoids an unnecessary low-precision intermediate during the conversion; it does not invent additional original detail. Black maps to video black (Y=16, Cb=Cr=128). No HDR metadata is added.

For comparison with original sRGB PNGs, decode with `BT709_TO_SRGB_FILTER` from the helper. A simple YCbCr-to-RGB matrix conversion alone does not invert the display transfer. Run `verify_instagram_export.py` to reproduce the synthetic patch audit recorded in `instagram-export-audit.json`; final file/loop checks are separate. These checks assess our export, not the quality of an actual Instagram upload.

## What the linked discussion supports

The [Reddit thread supplied for this project](https://www.reddit.com/r/editors/comments/16tunxo/instagram_crunched_up_my_videos_quality_how_do_i/) recommends enabling **Upload at highest quality**, and includes a report of good results exporting a 4K source into a 1080 × 1920 sequence at around 20 Mbps. Other commenters recommend below 10 Mbps, around 10 Mbps, or HEVC at 30 Mbps. Those reports conflict and are personal observations, not controlled comparisons or official requirements. The discussion also correctly distinguishes high-resolution source capture/supersampling from the resolution of the final delivery.

We do not rely on the thread's claim that small files avoid Instagram recompression. [Meta's March 2026 engineering description](https://engineering.fb.com/2026/03/02/video-engineering/ffmpeg-at-meta-media-processing-at-scale/) explains that uploaded videos generate multiple delivery encodings and that playback can select versions with different resolutions, codecs, frame rates and quality levels according to conditions. No local export setting guarantees a particular delivered version.

[Meta's official Instagram Reels API sample](https://github.com/fbsamples/reels_publishing_apis/blob/main/insta_reels_publishing_api_sample/README.md) documents MP4/MOV, H.264 or HEVC, progressive scan, closed GOP, 4:2:0, 23–60 fps, a recommended 9:16 ratio, a maximum horizontal size of 1920 pixels, and a maximum variable bitrate of 25 Mbps. These are API publishing specifications, not a claim that 25 Mbps is an optimal upload bitrate. The chosen delivery fits within them. [Meta's color-pipeline explanation](https://engineering.fb.com/2023/07/17/video-engineering/hdr-video-reels-meta/) describes HDR/SDR tone mapping and adaptive delivery; our original SDR content is kept SDR.

## Posting the file

1. Transfer the original MP4 to your phone using a method that preserves the file, such as AirDrop or a file download. Avoid an extra export or a messaging app's compressed-video mode
2. In Instagram's settings, find **Media quality** and enable **Upload at highest quality** if that control is available; exact menu wording can vary by app version
3. Select the original file, keep its vertical framing, and check the placement preview for the title and account name before posting
4. Copy the supplied `.text` description, including its model assumptions, study credits, and code link, into the Instagram description
5. Check the title and account name against the interface overlays. [Meta's Reels guidance](https://www.facebook.com/business/ads/facebook-instagram-reels-ads) recommends safe zones for key messages, but organic Reel controls, expanded captions, and different screens can cover different areas. The full-frame simulation continues beneath those areas; no lower blank panel is reserved

The upload-quality setting comes from the supplied discussion; the Instagram Help page and detailed Ads Guide were not retrievable during this review. No universal current pixel-safe area or guaranteed post-upload sharpness is claimed. The final appearance still depends on Instagram's processing, the viewer's device and playback conditions.

Technical reference: [FFmpeg zscale color-conversion documentation](https://ffmpeg.org/ffmpeg-filters.html#zscale-1)

Transfer reference: [zimg's display/scene transfer-function implementation](https://github.com/sekrit-twc/zimg/blob/master/src/zimg/colorspace/gamma.cpp)
