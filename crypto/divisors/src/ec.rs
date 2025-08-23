//! Trait for efficient coordinate conersion and general implementation
//! based on projective arithmetic.

use core::ops::Neg;
use std_shims::vec::Vec;

use zeroize::Zeroize;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};
use ff::Field;

use crate::DivisorCurve;

/// Point for which batch conversion to Weierstrass (X,Y) is cheap.
/// May and should be trivial for most curves.
/// Curve params are provided in case they are not available as constants.
pub trait XyPoint<F: Field>:
  Sized + Clone + Copy + Neg<Output = Self> + Zeroize + ConstantTimeEq + ConditionallySelectable
{
  /// The identity.
  const IDENTITY: Self;
  /// Add.
  fn add(x1: Self, x2: Self, curve: &Curve<F>) -> Self;
  /// Double.
  fn double(self, curve: &Curve<F>) -> Self;
  /// Efficient batch to_xy conversion.
  fn to_xy_batched(points: Vec<Self>) -> Vec<(F, F)>;
  /// If this point is the identity.
  fn is_identity(&self) -> Choice;
}

/// A point in projective coordinates
#[derive(Debug, Clone, Copy)]
pub struct Projective<F: Field> {
  x: F,
  y: F,
  z: F,
}

impl<F: Field> Projective<F> {
  /// (x,y,z)
  fn coordinates(self) -> (F, F, F) {
    let Self { x, y, z } = self;
    (x, y, z)
  }
  #[cfg(test)]
  fn to_affine_slow(self) -> (F, F) {
    let Self { x, y, z } = self;
    let z_inv = z.invert().unwrap();
    (x * z_inv, y * z_inv)
  }
}

impl<F: Field + Zeroize> Zeroize for Projective<F> {
  fn zeroize(&mut self) {
    let Self { x, y, z } = self;
    x.zeroize();
    y.zeroize();
    z.zeroize();
  }
}

impl<F: Field> ConstantTimeEq for Projective<F> {
  fn ct_eq(&self, other: &Self) -> Choice {
    let c1 = (self.x * other.z).ct_eq(&(other.x * self.z));
    let c2 = (self.y * other.z).ct_eq(&(other.y * self.z));
    c1 & c2
  }
}

impl<F: Field> Neg for Projective<F> {
  type Output = Self;
  fn neg(self) -> Self::Output {
    let Self { x, y, z } = self;
    let y = -y;
    Self { x, y, z }
  }
}

impl<F: Field> ConditionallySelectable for Projective<F> {
  fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
    let x = F::conditional_select(&a.x, &b.x, choice);
    let y = F::conditional_select(&a.y, &b.y, choice);
    let z = F::conditional_select(&a.z, &b.z, choice);
    Projective { x, y, z }
  }
}

impl<C: DivisorCurve> From<C> for Projective<C::FieldElement> {
  fn from(point: C) -> Self {
    let Some((x, y)) = C::to_xy(point) else { return Self::IDENTITY };
    Self { x, y, z: C::FieldElement::ONE }
  }
}

/// Curve params x^3 + ax + b
pub struct Curve<F> {
  /// A in x^3 + Ax + B
  pub a: F,
  /// B in x^3 + Ax + B
  pub b: F,
}

impl<F: Field + Zeroize> XyPoint<F> for Projective<F> {
  const IDENTITY: Self = Projective { x: F::ZERO, y: F::ONE, z: F::ZERO };

  // based on 13.2.1.b of https://hyperelliptic.org/HEHCC
  fn add(x1: Self, x2: Self, curve: &Curve<F>) -> Self {
    let identity_left = CtOption::new(x1, x2.is_identity());
    let identity_right = CtOption::new(x2, x1.is_identity());
    let double = x1.double(curve);
    let equals = x1.ct_eq(&x2);
    let (x1, y1, z1) = x1.coordinates();
    let (x2, y2, z2) = x2.coordinates();
    let a = y2 * z1 - y1 * z2;
    let b = x2 * z1 - x1 * z2;
    let z1z2 = z1 * z2;
    let bb = b.square();
    let bbb = bb * b;
    let c = a.square() * z1z2 - bbb - bb.double() * x1 * z2;

    let x = b * c;
    let y = a * (bb * x1 * z2 - c) - bbb * y1 * z2;
    let z = bbb * z1z2;
    let add = Projective { x, y, z };
    let res = Self::conditional_select(&add, &double, equals);
    identity_left.unwrap_or(identity_right.unwrap_or(res))
  }

  fn double(self, curve: &Curve<F>) -> Self {
    let (x1, y1, z1) = self.coordinates();
    let xx = x1.square();
    let xx3 = xx.double() + xx;
    let a = curve.a * z1.square() + xx3;
    let b = y1 * z1;
    let c = x1 * y1 * b;
    let c4 = c.double().double();
    let c8 = c4.double();
    let d = a.square() - c8;
    let x = (b * d).double();
    let bb = b.square();
    let bb8 = bb.double().double().double();
    let y = a * (c4 - d) - bb8 * y1.square();
    let z = bb8 * b;
    Self { x, y, z }
  }

  fn to_xy_batched(points: Vec<Self>) -> Vec<(F, F)> {
    let mut z = points.iter().map(|p| p.z).collect::<Vec<F>>();
    let mut scratch_space = vec![F::ZERO; z.len()];
    ff::BatchInverter::invert_with_external_scratch(&mut z, &mut scratch_space);
    points.into_iter().zip(z).map(|(p, z_inv)| (p.x * z_inv, p.y * z_inv)).collect()
  }

  fn is_identity(&self) -> Choice {
    self.z.ct_eq(&F::ZERO)
  }
}

#[cfg(test)]
mod ed25519_test {
  use group::Group;
  use dalek_ff_group::EdwardsPoint;
  use crate::{
    ec::{Curve, Projective, XyPoint},
    DivisorCurve,
  };

  #[test]
  fn projective() {
    let a = EdwardsPoint::a();
    let b = EdwardsPoint::b();
    let curve = Curve { a, b };
    let to_xy = |p| EdwardsPoint::to_xy(p).unwrap();

    let point = EdwardsPoint::generator();
    let projective = Projective::from(point);
    assert_eq!(projective.to_affine_slow(), to_xy(point));
    // println!("P = {:?}", projective.to_affine_slow());

    let doubled = point.double();
    // println!("2P = {:?}", to_xy(doubled));
    let doubled_projective = projective.double(&curve);
    // println!("2P = {:?}", doubled_projective.to_affine_slow());
    assert_eq!(doubled_projective.to_affine_slow(), to_xy(doubled));

    let triple = doubled + point;
    // println!("3P = {:?}", to_xy(triple));
    let triple_projective = Projective::add(doubled_projective, projective, &curve);
    // println!("3P = {:?}", triple_projective.to_affine_slow());
    assert_eq!(triple_projective.to_affine_slow(), to_xy(triple));
  }
}
