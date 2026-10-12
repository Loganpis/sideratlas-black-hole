# Spectral color data

The unchanged `CIE_xyz_1931_2deg.csv` and its metadata were retrieved from the
International Commission on Illumination on 2026-10-07.

CIE 2019, *Colour-matching functions of CIE 1931 standard colorimetric observer*,
International Commission on Illumination (CIE), Vienna, AT.
DOI: [10.25039/CIE.DS.xvudnb9b](https://doi.org/10.25039/CIE.DS.xvudnb9b).

- [Dataset description](https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer)
- [Official CSV](https://files.cie.co.at/Publications-datasets/CIE_xyz_1931_2deg.csv)
- [Official metadata](https://files.cie.co.at/Publications-datasets/CIE_xyz_1931_2deg.csv_metadata.json)
- Official MD5: `17cca777db64b17170f06f67ce9d3ab7`

Source publication: CIE 018:2019, Table 6; 360–830 nm inclusive at 1 nm steps.
License and provenance are preserved in the supplied official metadata file.

The integration computes the Planck spectrum at the observed temperature
`T_obs = g T_eff`, integrates against all three CIE matching functions, and
transforms XYZ to linear sRGB using section 7 of the
[ICC sRGB interpretation](https://registry.color.org/rgb-registry/files/sRGB.pdf).
There is no arbitrary blackbody-color approximation and no normalization of
individual pixels. Display exposure and dynamic-range compression are distinct
operations in the renderer. The integration retains relative visible brightness;
the common photometric scale of 683 lm/W is omitted because it is absorbed by
the one fixed display exposure.
