use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Rem, Sub, SubAssign},
};

use num::{Bounded, One, ToPrimitive, Zero};

pub trait Data = Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + MulAssign
    + Div<Output = Self>
    + DivAssign
    + Rem
    + One
    + Zero
    + PartialEq
    + PartialOrd
    + Bounded
    + Copy
    + Clone
    + Debug
    + Default
    + ToPrimitive
    + Send
    + Unpin
    + 'static;
