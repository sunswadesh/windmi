# windmi
A low order solar wind magnetosphere ionosphere model. Modular machine learning tools are used for parameter optimization, but in principal any algorithm can be used.

**Model Developer(s)**
W. Horton, M. L. Mays, E. Spencer and I. Doxas
Institute for Fusion Studies, Department of Physics, University of Texas at Austin

**Maintained and modified by:**
Swadesh Patra

**Versions: Model adaptations:**
1. Magnetospheric currents combined to estimate low-latitude ground magnetic perturbations (Predicting SymH index).
2. Variable ionospheric conductivity.
3. Solar wind dependent ring current decay time
4. Tail current validation. 



Original Model Description
WINDMI is a low-dimensional model of the energy transfer from the solar wind through the magnetosphere and into the ionosphere. The model uses the analogy of electric circuitry (capacitances, resistances, inductances) to describe with a set of nonlinear ordinary differential equations the response of the magnetosphere-ionosphere system to solar wind driving. The electric driving voltage applied by the solar wind is described either by the Rectified (Reiff and Luhmann, 1986), Siscoe (Siscoe et al. 2002) or Newell (2007) coupling function (references in Spencer et al. 2009).

Model Input
The model can be driven by ACE Real-Time or Level2 solar wind measurements over a few days and provides as outputs various magnetospheric and ionospheric parameters (region I and II currents, cross-polarcap voltages, ionospheric dissipation).

Model Output
The major products of the model are the energy of the ring current (Dst/SymH index) and the auroral indices(AU,AL) of geomagnetic activity.

References and relevant publications
Horton, W., and I. Doxas (1998), A low-dimensional dynamical model for the  solar wind  driven geotail-ionosphere  system, J.  Geophys. Res., 103(A3), 4561-4572.

Spencer, E.,  W. Horton, M. L.  Mays, I. Doxas, and  J. Kozyra (2007), Analysis  of the  3-7 October  2000 and  15-24 April  2002 geomagnetic storms with an optimized  nonlinear dynamical model, J. Geophys. Res., 112, A04S90, doi:10.1029/2006JA012019.
