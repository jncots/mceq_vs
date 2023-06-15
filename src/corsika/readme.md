This is a message from Tetiana Kozynets about initial structure
of the file with Corsika results. The file format is the base for later changes and imporvments:

```
is it OK if I give you this file? 
Anton can then bin the energies :slightly_smiling_face:
the keys of the HDF5 file are PDG IDs (±12, ±13, ±14), 
as well as 'num_primaries' which tells you the number 
of simulated protons. for each PDG ID key there are 4 
observation levels (1, 2, 3, 4), which correspond to 
143, 638, 1167, and 1195 g/cm2. then for each PDG key 
and observation level there are two columns — energy and theta.
so for example, if you load this file under the name 
corsika_data, all muon energies at 143 g/cm2 can be 
extracted as np.array(corsika_data['13']['1']['energy [GeV]']) . 
and then num_primaries can be used to normalize it
```