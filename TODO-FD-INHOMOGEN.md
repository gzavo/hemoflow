# Flow diverter: device-level, area-weighted resistance (plan)

Status: planned (November 2026). Per-voxel resistance multipliers are in place
(TODO.md, section 2). This plan changes what the coefficients mean, from one voxel
layer to the whole device, so that a spatially varying resistance is no longer
dominated by how the device happens to be voxelized.

## Why the current model is not enough

In every porous voxel the solver applies the force density

    f = -(m_lin a + m_quad b |u|) u / dx

so `linCoeff` (a) and `quadCoeff` (b) are the pressure drop across **one voxel layer**,
and the resistance of the device is proportional to the number of layers a fluid
particle crosses.

- **That number depends on the voxelization.** Device voxels counted along the device
  normal through each triangle centre (±3 dx), weighted by triangle area, for the
  NAP180 + PED 5×20 case (current preprocessor):

  | dx | 0 voxels (gap) | 1 voxel | 2 voxels | 3 voxels |
  |---|---|---|---|---|
  | 0.36 mm (1 M voxels) | 6.5 % | 78.3 % | 14.5 % | 0.7 % |
  | 0.21 mm (5 M voxels) | 7.2 % | 76.1 % | 16.2 % | 0.4 % |
  | 0.13 mm (20 M voxels) | 3.0 % | 77.1 % | 18.5 % | 1.4 % |

  The mix is about the same at every dx, so the resistance of the device as a whole
  hardly depends on dx. Locally it does depend on how the surface is oriented to the
  grid axes: where the device is crossed as two voxels (15–19 % of its area) the
  resistance is doubled by the voxelization alone, in patches that move with the
  device position and the grid. The `Density` descriptor of the shipped device varies
  0.83–1.78 over 90 % of its points (0.57–2.26 in total), a variation of the same size.
  Geometry files from the preprocessor before November 2026 (shifted projections)
  had 50–52 % at one voxel and 35–40 % at two.
- **Gaps.** On 3–7 % of the device area no device voxel lies along the normal, so flow
  crosses there without resistance: the axis-aligned rays miss parts of the surface.
- **Interpolation.** Each device voxel takes the values of the nearest of the 812 points
  of `<stent_mesh_base>values.vtp` (about one per braid cell, ~0.5 mm), which gives a
  patchwork at the voxel scale.

## Target model

- a and b describe the **whole device**: Δp = a u + b u² across it, per unit device area.
  This is what a bench test measures, and it does not depend on dx.
- Each voxel gets a weight w = A_v / dx², where A_v is the device surface area inside
  the voxel. The force density becomes

      f = -w (m_lin a u + m_quad b |u| u) / dx

  so the force summed over the voxels of a device patch of area A is
  Σ f dx³ = (a u + b |u| u) Σ A_v = (a u + b |u| u) A, whatever the orientation of the
  device to the grid and whatever dx is.
- The quadratic term stays b |u| u (rotation invariant since November 2026).
- Optionally the resistance acts mainly on the flow **through** the device: a per-voxel
  normal n (area-weighted mean of the triangle normals in the voxel), the full
  resistance for the normal velocity u_n and a factor κ_t (0 ≤ κ_t ≤ 1) for the
  tangential velocity. Physics decision, see below.

## Preprocessor

- [ ] Subdivide the device triangles until every edge is shorter than dx/4. Each
      sub-triangle carries its area, its centroid and the descriptors interpolated
      barycentrically from the vertices of its parent triangle.
- [ ] Bin the sub-triangles into voxels by centroid: per voxel the area A_v (and, if
      anisotropic, the area-weighted normal), the area-weighted mean descriptors.
      Numpy `bincount` on flat voxel indices; at dx = 0.13 mm the PED mesh gives about
      100 k sub-triangles, so this costs well under a second.
- [ ] Keep only lumen voxels, as now (area in wall and opening voxels is dropped; the
      device is apposed to the wall there). Print the binned area against the STL area
      inside the lumen as a check.
- [ ] Every voxel that contains device area becomes porous, so there are no gaps: a
      lattice link (also a diagonal one) between two cell centres lies in the union of
      its two end voxels, so any link that crosses the surface has a porous end. This
      also replaces the ray-cast shell voxelization (TODO.md section 2, item 2).
- [ ] Map the descriptors to multipliers (function in `main.py`, see "Physics inputs").
- [ ] Write new arrays, so that a file says which model it uses: `deviceWeightLinear`
      and `deviceWeightQuadratic` = w · m_lin and w · m_quad (float32), and optionally
      `deviceNormal` (float32, (Nx, Ny, Nz, 3)). Add them to `geometry.vti`.

## Solver

- [ ] Read the new arrays like the multipliers (geometry.cpp, `scatterComponents`).
      Two components as now, five with the normal.
- [ ] `PorousForceFunctional`: f = -(W_lin a u + W_quad b |u| u) / dx with the per-voxel
      weights W, plus the normal/tangential split if anisotropic.
- [ ] The geometry file decides the model: `stent` + multipliers = layer model (as
      now), the new arrays = device model; both in one file is an error. The log states
      the model and what `linCoeff` and `quadCoeff` mean.
- [ ] Optional: allocate the weight field only over the bounding box of the device
      instead of the whole lattice.
- [ ] Optional: an output field with the porous force density, to check the momentum
      sink directly.

## Physics inputs needed

1. What `Density`, `Expansion` and `Force` in `PED_5x20_values.vtp` mean (relative to
   the nominal, unconstrained device?).
2. a and b (or porosity and permeability) as functions of these descriptors: a formula
   or a fit to bench data.
3. Isotropic resistance, or mainly through the device (κ_t).
4. Carreau runs: whether the Darcy term should follow the local viscosity
   (a = μ t / K), which needs a permeability K rather than a fixed a.

## Verification

- [ ] **Momentum sink vs surface integral.** A flat porous disk across the straight
      pipe of `tests/verification/0-pipe`, tilted by 0°, 20°, 35° and 45°, at three dx,
      Darcy term only. The total porous force from the solution must equal
      a ∫ u dA over the device STL (velocity interpolated on the surface) to within a
      few percent for every angle and dx. The layer model run alongside shows its
      spread.
- [ ] **Pressure drop.** For the perpendicular disk, the jump in cross-section mean
      pressure (with the Poiseuille gradient removed) equals a Q/A + b ∫u² dA / A.
- [ ] **Two regions.** Disk halves with multipliers 1 and 4: the flux split and the
      pressure drop follow from Δp = m a u per half.
- [ ] **Grid convergence on NAP180** (1 M, 5 M, 20 M voxels): the inflow into the
      aneurysm should converge with the device model; record the layer model alongside.
- [ ] **Device position.** Shift the device STL by fractions of a voxel (0, 1/4, 1/2 dx):
      the aneurysm inflow must not change with the device model, while the patches of
      double resistance move with the layer model.
- [ ] Against bench measurements, if available.

## Compatibility and cost

- Older geometry files and existing configs keep the layer model. Coefficients tuned
  for the layer model (e.g. 100 / 200 in `testcase/`) do not carry over; the device
  model needs device-level values.
- Solver cost as now (two components); the normal adds three components per cell.
- Stability: the explicit porous force is fine while w m a dt / (ρ dx) stays well
  below 1. With the shipped coefficients it is about 0.005 at 5 M voxels (dx 0.21 mm,
  dt 1e-5 s); a semi-implicit force (TODO.md section 2) is only needed for much more
  resistive devices.

## How the table above was measured

For every triangle of the device STL whose centre lies in a lumen voxel: sample the
line along the triangle normal through its centre over ±3 dx, collect the voxels it
passes through and count those marked in `stent`; weight by the triangle area.
