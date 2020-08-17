# HemoFlow

Macroscopic flow simulation, aimed at vessel simulations. It relies on the preprocess_LBM package to prepare the simulation domain.

# Note
- All input parameters should be SI or non-dimensional!

## Shortcomings
- Openings must be on the axis aligned (AA) bounding box border for now to make geometry preparation automatic
- An opening cannot fall to an edge or corner of the AA bounding box


## TODO
- New inlet scale function (a more realistic one)
- New outlet pressure distribution based on Murray-law (smallest outlet -> p=0)
- Pass the angle of the openings based on centerline calculations (new ID / opening, every centerline goes from the inlet to an opening)

