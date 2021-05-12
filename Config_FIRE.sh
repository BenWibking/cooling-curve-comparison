#!/bin/bash            # this line only there to enable syntax highlighting in this file

####################################################################################################
# -------------------------------------- Radiative Cooling physics (mostly geared towards galactic/extragalactic cooling)
# -------------------------- These modules were originally developed for a combination of proprietary physics modules. However they are now written in
# --------------------------   a form which allows them to be modular (and public). Users are free to use the Grackle modules and standard 'COOLING' flags,
# --------------------------   provided proper credit/citations are provided to the relevant methods papers given in the Users Guide ---
# --------------------------   but all users should cite Hopkins et al. 2017 (arXiv:1702.06148), where Appendix B details the cooling physics
OUTPUT_COOLRATE_DETAIL
METALS                         # enable metallicities (with multiple species optional) for gas and stars [must be included in ICs or injected via dynamical feedback; needed for some routines]
COOLING                        # enables radiative cooling and heating: if GALSF, also external UV background read from file "TREECOOL" (included in the cooling folder; be sure to cite its source as well, given in the TREECOOL file)
COOL_LOW_TEMPERATURES          # allow fine-structure and molecular cooling to ~10 K; account for optical thickness and line-trapping effects with proper opacities [requires METALS]
COOL_LOWTEMP_THIN_ONLY
COOL_METAL_LINES_BY_SPECIES    # use full multi-species-dependent cooling tables ( http://www.tapir.caltech.edu/~phopkins/public/spcool_tables.tgz, or the Bitbucket site); requires METALS on; cite Wiersma et al. 2009 (MNRAS, 393, 99) in addition to Hopkins et al. 2017 (arXiv:1702.06148)
#COOL_GRACKLE                   # enable Grackle: cooling+chemistry package (requires COOLING above; https://grackle.readthedocs.org/en/latest ); see Grackle code for their required citations
#COOL_GRACKLE_CHEMISTRY=1       # choose Grackle cooling chemistry: (0)=tabular, (1)=Atomic, (2)=(1)+H2+H2I+H2II, (3)=(2)+DI+DII+HD
## ----------------------------------------------------------------------------------------------------

####################################################################################################
# ----------------- Galaxy formation & Galactic Star formation
####################################################################################################
## ---------------------------------------------------------------------------------------------------
GALSF                           # master switch for galactic star formation model: enables SF, stellar ages, generations, etc. [cite Springel+Hernquist 2003, MNRAS, 339, 289]
