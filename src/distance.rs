//! Distances on the Earth.
//!
//! * [`haversine`]: great-circle distance on a **sphere** of radius
//!   [`EARTH_MEAN_RADIUS_M`]. Closed-form, fast, always defined. Because the
//!   Earth is an oblate ellipsoid, the result differs from the true geodesic
//!   by up to roughly 0.5% (about 0.56% for short north-south hops near the
//!   equator, where the meridional radius of curvature is smallest).
//! * [`vincenty_inverse`]: Vincenty's (1975) iterative solution of the inverse
//!   geodesic problem on the **WGS84 ellipsoid**. Accurate to well under a
//!   millimetre when it converges, but the iteration on `lambda` can fail to
//!   converge for nearly antipodal points; that case is reported as
//!   [`VincentyError::NoConvergence`] instead of returning a wrong number.

use std::f64::consts::PI;
use std::fmt;

use crate::types::LonLat;

/// IUGG mean Earth radius `R1 = (2a + b) / 3`, in metres.
pub const EARTH_MEAN_RADIUS_M: f64 = 6_371_008.8;

/// WGS84 ellipsoid parameters.
pub mod wgs84 {
    /// Semi-major axis (equatorial radius), metres.
    pub const A: f64 = 6_378_137.0;
    /// Flattening.
    pub const F: f64 = 1.0 / 298.257_223_563;
    /// Semi-minor axis (polar radius), metres.
    pub const B: f64 = A * (1.0 - F);
}

/// Maximum number of `lambda` iterations before giving up.
pub const VINCENTY_MAX_ITERATIONS: u32 = 200;
/// Convergence threshold on `lambda` (radians). 1e-12 rad is ~0.006 mm.
pub const VINCENTY_TOLERANCE: f64 = 1e-12;

/// Great-circle distance in metres using the Haversine formula on a sphere of
/// mean Earth radius.
pub fn haversine(a: LonLat, b: LonLat) -> f64 {
    haversine_with_radius(a, b, EARTH_MEAN_RADIUS_M)
}

/// Haversine distance on a sphere of the given radius (same unit as the result).
pub fn haversine_with_radius(a: LonLat, b: LonLat, radius: f64) -> f64 {
    let phi1 = a.lat.to_radians();
    let phi2 = b.lat.to_radians();
    let d_phi = (b.lat - a.lat).to_radians();
    let d_lambda = (b.lon - a.lon).to_radians();

    let s_phi = (d_phi * 0.5).sin();
    let s_lambda = (d_lambda * 0.5).sin();
    // h = hav(central angle); clamp guards against rounding slightly above 1
    // for antipodal points, which would make asin return NaN.
    let h = (s_phi * s_phi + phi1.cos() * phi2.cos() * s_lambda * s_lambda).clamp(0.0, 1.0);
    2.0 * radius * h.sqrt().asin()
}

/// Initial bearing (forward azimuth) on the sphere, degrees in `[0, 360)`.
pub fn initial_bearing_sphere(a: LonLat, b: LonLat) -> f64 {
    let phi1 = a.lat.to_radians();
    let phi2 = b.lat.to_radians();
    let d_lambda = (b.lon - a.lon).to_radians();
    let y = d_lambda.sin() * phi2.cos();
    let x = phi1.cos() * phi2.sin() - phi1.sin() * phi2.cos() * d_lambda.cos();
    normalize_bearing(y.atan2(x).to_degrees())
}

fn normalize_bearing(deg: f64) -> f64 {
    let b = deg.rem_euclid(360.0);
    // rem_euclid can return exactly 360.0 for tiny negative inputs.
    if b >= 360.0 {
        0.0
    } else {
        b
    }
}

/// Result of a successful Vincenty inverse computation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Geodesic {
    /// Ellipsoidal distance in metres.
    pub distance_m: f64,
    /// Forward azimuth at the start point, degrees clockwise from north, `[0, 360)`.
    pub initial_bearing_deg: f64,
    /// Forward azimuth at the end point (direction of travel on arrival), `[0, 360)`.
    pub final_bearing_deg: f64,
    /// Number of `lambda` iterations used.
    pub iterations: u32,
}

/// Why [`vincenty_inverse`] could not produce a result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VincentyError {
    /// The iteration did not converge (typically nearly antipodal points).
    NoConvergence { iterations: u32 },
}

impl fmt::Display for VincentyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VincentyError::NoConvergence { iterations } => write!(
                f,
                "Vincenty inverse did not converge after {iterations} iterations (points nearly antipodal?)"
            ),
        }
    }
}

impl std::error::Error for VincentyError {}

/// Vincenty inverse on WGS84: distance and azimuths between two positions.
pub fn vincenty_inverse(p1: LonLat, p2: LonLat) -> Result<Geodesic, VincentyError> {
    vincenty_inverse_on(p1, p2, wgs84::A, wgs84::F)
}

/// Convenience wrapper returning only the distance in metres.
pub fn vincenty_distance(p1: LonLat, p2: LonLat) -> Result<f64, VincentyError> {
    vincenty_inverse(p1, p2).map(|g| g.distance_m)
}

/// Vincenty inverse on an arbitrary ellipsoid (`a` semi-major axis, `f` flattening).
pub fn vincenty_inverse_on(
    p1: LonLat,
    p2: LonLat,
    a: f64,
    f: f64,
) -> Result<Geodesic, VincentyError> {
    let b = a * (1.0 - f);

    // Longitude difference, normalised into [-pi, pi].
    let l = {
        let d = (p2.lon - p1.lon).to_radians();
        (d + PI).rem_euclid(2.0 * PI) - PI
    };

    // Reduced latitudes (latitude on the auxiliary sphere).
    let u1 = ((1.0 - f) * p1.lat.to_radians().tan()).atan();
    let u2 = ((1.0 - f) * p2.lat.to_radians().tan()).atan();
    let (sin_u1, cos_u1) = u1.sin_cos();
    let (sin_u2, cos_u2) = u2.sin_cos();

    let mut lambda = l;
    for iteration in 1..=VINCENTY_MAX_ITERATIONS {
        let (sin_lambda, cos_lambda) = lambda.sin_cos();
        let t1 = cos_u2 * sin_lambda;
        let t2 = cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda;
        let sin_sigma = (t1 * t1 + t2 * t2).sqrt();
        let cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda;

        if sin_sigma == 0.0 {
            if cos_sigma > 0.0 {
                // Coincident points.
                return Ok(Geodesic {
                    distance_m: 0.0,
                    initial_bearing_deg: 0.0,
                    final_bearing_deg: 0.0,
                    iterations: iteration,
                });
            }
            // Exactly antipodal on the auxiliary sphere: azimuth undefined.
            return Err(VincentyError::NoConvergence {
                iterations: iteration,
            });
        }

        let sigma = sin_sigma.atan2(cos_sigma);
        let sin_alpha = cos_u1 * cos_u2 * sin_lambda / sin_sigma;
        let cos_sq_alpha = 1.0 - sin_alpha * sin_alpha;
        // On the equator cos^2(alpha) = 0 and cos(2 sigma_m) is irrelevant.
        let cos_2sigma_m = if cos_sq_alpha != 0.0 {
            cos_sigma - 2.0 * sin_u1 * sin_u2 / cos_sq_alpha
        } else {
            0.0
        };
        let c = f / 16.0 * cos_sq_alpha * (4.0 + f * (4.0 - 3.0 * cos_sq_alpha));
        let lambda_prev = lambda;
        lambda = l
            + (1.0 - c)
                * f
                * sin_alpha
                * (sigma
                    + c * sin_sigma
                        * (cos_2sigma_m + c * cos_sigma * (-1.0 + 2.0 * cos_2sigma_m * cos_2sigma_m)));

        if !lambda.is_finite() {
            return Err(VincentyError::NoConvergence {
                iterations: iteration,
            });
        }

        if (lambda - lambda_prev).abs() < VINCENTY_TOLERANCE {
            let u_sq = cos_sq_alpha * (a * a - b * b) / (b * b);
            let big_a =
                1.0 + u_sq / 16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
            let big_b = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));
            let c2m2 = cos_2sigma_m * cos_2sigma_m;
            let delta_sigma = big_b
                * sin_sigma
                * (cos_2sigma_m
                    + big_b / 4.0
                        * (cos_sigma * (-1.0 + 2.0 * c2m2)
                            - big_b / 6.0
                                * cos_2sigma_m
                                * (-3.0 + 4.0 * sin_sigma * sin_sigma)
                                * (-3.0 + 4.0 * c2m2)));
            let distance_m = b * big_a * (sigma - delta_sigma);

            let (sin_lambda, cos_lambda) = lambda.sin_cos();
            let alpha1 = (cos_u2 * sin_lambda).atan2(cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda);
            let alpha2 =
                (cos_u1 * sin_lambda).atan2(-sin_u1 * cos_u2 + cos_u1 * sin_u2 * cos_lambda);

            return Ok(Geodesic {
                distance_m,
                initial_bearing_deg: normalize_bearing(alpha1.to_degrees()),
                final_bearing_deg: normalize_bearing(alpha2.to_degrees()),
                iterations: iteration,
            });
        }
    }

    Err(VincentyError::NoConvergence {
        iterations: VINCENTY_MAX_ITERATIONS,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dms(d: f64, m: f64, s: f64) -> f64 {
        LonLat::dms_to_deg(d, m, s)
    }

    /// Test case from Vincenty (1975), also in Geoscience Australia's docs:
    /// Flinders Peak -> Buninyong.
    fn flinders_buninyong() -> (LonLat, LonLat) {
        let flinders = LonLat::new(dms(144.0, 25.0, 29.52440), dms(-37.0, 57.0, 3.72030));
        let buninyong = LonLat::new(dms(143.0, 55.0, 35.38390), dms(-37.0, 39.0, 10.15610));
        (flinders, buninyong)
    }

    const ALGIERS: LonLat = LonLat::new(3.0588, 36.7538);
    const ORAN: LonLat = LonLat::new(-0.6308, 35.6971);

    #[test]
    fn vincenty_flinders_peak_buninyong() {
        let (a, b) = flinders_buninyong();
        let g = vincenty_inverse(a, b).unwrap();
        assert!((g.distance_m - 54_972.271).abs() < 0.001, "{}", g.distance_m);
        // Forward azimuth 306° 52' 05.37"
        assert!((g.initial_bearing_deg - dms(306.0, 52.0, 5.37)).abs() < 1.0 / 3600.0 * 0.01);
        // Reverse azimuth 127° 10' 25.07" => forward azimuth at end = 307° 10' 25.07"
        assert!((g.final_bearing_deg - dms(307.0, 10.0, 25.07)).abs() < 1.0 / 3600.0 * 0.01);
        assert!(g.iterations < 10);
    }

    #[test]
    fn vincenty_equator_one_degree() {
        // Along the equator the geodesic is the equator: a * (pi / 180).
        let g = vincenty_inverse(LonLat::new(0.0, 0.0), LonLat::new(1.0, 0.0)).unwrap();
        assert!((g.distance_m - 111_319.490_793).abs() < 1e-5, "{}", g.distance_m);
        assert!((g.initial_bearing_deg - 90.0).abs() < 1e-9);
    }

    #[test]
    fn vincenty_quarter_meridian() {
        // Equator to pole along a meridian: 10 001 965.729 m on WGS84.
        let d = vincenty_distance(LonLat::new(0.0, 0.0), LonLat::new(0.0, 90.0)).unwrap();
        assert!((d - 10_001_965.729).abs() < 0.001, "{d}");
    }

    #[test]
    fn vincenty_coincident_points() {
        let g = vincenty_inverse(ALGIERS, ALGIERS).unwrap();
        assert_eq!(g.distance_m, 0.0);
    }

    #[test]
    fn vincenty_symmetry() {
        let d1 = vincenty_distance(ALGIERS, ORAN).unwrap();
        let d2 = vincenty_distance(ORAN, ALGIERS).unwrap();
        assert!((d1 - d2).abs() < 1e-6);
    }

    #[test]
    fn vincenty_antipodal_fails_cleanly() {
        // Exactly antipodal on the equator.
        let r = vincenty_inverse(LonLat::new(0.0, 0.0), LonLat::new(180.0, 0.0));
        assert!(matches!(r, Err(VincentyError::NoConvergence { .. })), "{r:?}");
        // Nearly antipodal: the classic failure region of the lambda iteration.
        let r = vincenty_inverse(LonLat::new(0.0, 0.0), LonLat::new(179.7, 0.5));
        assert!(matches!(r, Err(VincentyError::NoConvergence { .. })), "{r:?}");
        assert!(r.unwrap_err().to_string().contains("did not converge"));
    }

    #[test]
    fn algiers_oran() {
        let h = haversine(ALGIERS, ORAN);
        let v = vincenty_distance(ALGIERS, ORAN).unwrap();
        // About 352 km as the crow flies.
        assert!((345_000.0..360_000.0).contains(&h), "{h}");
        assert!((345_000.0..360_000.0).contains(&v), "{v}");
        assert!((h - v).abs() / v < 0.005);
        // Oran is west-south-west of Algiers.
        let bearing = vincenty_inverse(ALGIERS, ORAN).unwrap().initial_bearing_deg;
        assert!((240.0..260.0).contains(&bearing), "{bearing}");
    }

    #[test]
    fn haversine_known_values() {
        // One degree of arc on the mean sphere.
        let d = haversine(LonLat::new(0.0, 0.0), LonLat::new(0.0, 1.0));
        assert!((d - EARTH_MEAN_RADIUS_M * PI / 180.0).abs() < 1e-6);
        // Antipodes: half the circumference, no NaN.
        let d = haversine(LonLat::new(0.0, 0.0), LonLat::new(180.0, 0.0));
        assert!((d - EARTH_MEAN_RADIUS_M * PI).abs() < 1e-6);
        // Crossing the antimeridian is handled by the trigonometry.
        let d = haversine(LonLat::new(179.5, 0.0), LonLat::new(-179.5, 0.0));
        assert!((d - EARTH_MEAN_RADIUS_M * PI / 180.0).abs() < 1e-6);
        assert_eq!(haversine(ORAN, ORAN), 0.0);
    }

    #[test]
    fn haversine_vs_vincenty_flinders() {
        let (a, b) = flinders_buninyong();
        let h = haversine(a, b);
        assert!((h - 54_972.271).abs() / 54_972.271 < 0.005, "{h}");
    }

    #[test]
    fn vincenty_crosses_antimeridian() {
        let a = LonLat::new(179.5, 10.0);
        let b = LonLat::new(-179.5, 10.0);
        let d = vincenty_distance(a, b).unwrap();
        assert!(d < 120_000.0, "{d}");
    }

    #[test]
    fn sphere_bearing() {
        assert!((initial_bearing_sphere(LonLat::new(0.0, 0.0), LonLat::new(0.0, 10.0))).abs() < 1e-9);
        assert!((initial_bearing_sphere(LonLat::new(0.0, 0.0), LonLat::new(10.0, 0.0)) - 90.0).abs() < 1e-9);
        assert!((initial_bearing_sphere(LonLat::new(0.0, 0.0), LonLat::new(-10.0, 0.0)) - 270.0).abs() < 1e-9);
    }
}
