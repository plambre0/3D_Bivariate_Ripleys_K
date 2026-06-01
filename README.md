bivariate_k_3d(ref_points, eval_points, r_values, cube_bounds, precomputed_vols = NULL)//
Parameters\\
ref_points: Numeric matrix of reference (type-1) point coordinates. Each row is one point: columns are x, y, z. These are the 'from' points — spheres are centered here for counting.\\
eval_points: Numeric matrix of evaluation (type-2) point coordinates. Each row is one point: columns are x, y, z. These are the 'to' points that are counted inside each sphere.\\
r_values: Vector of search radii at which to evaluate K(r), L(r), and H(r). Values should be positive and typically span 0 to roughly half the shortest domain edge.\\
cube_bounds: Bounding box as c(xmin, xmax, ymin, ymax, zmin, zmax). Must enclose all points in both ref_points and eval_points.
precomputed_vols: Pre-computed clipped sphere volumes, one per entry of r_values, already summed over all ref_points. When NULL the volumes are computed internally via clipped_sphere_volume. Passing pre-computed values (e.g. from csr_envelope) avoids redundant calculation across simulation replicates.//
Returns\\
\\
Returns a vector of the used radii, Ripley's K, L, and H statistics.\\
