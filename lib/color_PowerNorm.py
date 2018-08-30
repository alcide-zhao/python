from collections import Sized
import itertools
import re
import warnings

import numpy as np
class Normalize(object):
    """
    A class which, when called, can normalize data into
    the ``[0.0, 1.0]`` interval.
    """
    def __init__(self, vmin=None, vmax=None, clip=False):
        """
        If *vmin* or *vmax* is not given, they are initialized from the
        minimum and maximum value respectively of the first input
        processed.  That is, *__call__(A)* calls *autoscale_None(A)*.
        If *clip* is *True* and the given value falls outside the range,
        the returned value will be 0 or 1, whichever is closer.
        Returns 0 if::
            vmin==vmax
        Works with scalars or arrays, including masked arrays.  If
        *clip* is *True*, masked values are set to 1; otherwise they
        remain masked.  Clipping silently defeats the purpose of setting
        the over, under, and masked colors in the colormap, so it is
        likely to lead to surprises; therefore the default is
        *clip* = *False*.
        """
        self.vmin = _sanitize_extrema(vmin)
        self.vmax = _sanitize_extrema(vmax)
        self.clip = clip

    @staticmethod
    def process_value(value):
        """
        Homogenize the input *value* for easy and efficient normalization.
        *value* can be a scalar or sequence.
        Returns *result*, *is_scalar*, where *result* is a
        masked array matching *value*.  Float dtypes are preserved;
        integer types with two bytes or smaller are converted to
        np.float32, and larger types are converted to np.float64.
        Preserving float32 when possible, and using in-place operations,
        can greatly improve speed for large arrays.
        Experimental; we may want to add an option to force the
        use of float32.
        """
        is_scalar = not cbook.iterable(value)
        if is_scalar:
            value = [value]
        dtype = np.min_scalar_type(value)
        if np.issubdtype(dtype, np.integer) or dtype.type is np.bool_:
            # bool_/int8/int16 -> float32; int32/int64 -> float64
            dtype = np.promote_types(dtype, np.float32)
        # ensure data passed in as an ndarray subclass are interpreted as
        # an ndarray. See issue #6622.
        mask = np.ma.getmask(value)
        data = np.asarray(np.ma.getdata(value))
        result = np.ma.array(data, mask=mask, dtype=dtype, copy=True)
        return result, is_scalar

    def __call__(self, value, clip=None):
        """
        Normalize *value* data in the ``[vmin, vmax]`` interval into
        the ``[0.0, 1.0]`` interval and return it.  *clip* defaults
        to *self.clip* (which defaults to *False*).  If not already
        initialized, *vmin* and *vmax* are initialized using
        *autoscale_None(value)*.
        """
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        # Convert at least to float, without losing precision.
        (vmin,), _ = self.process_value(self.vmin)
        (vmax,), _ = self.process_value(self.vmax)
        if vmin == vmax:
            result.fill(0)   # Or should it be all masked?  Or 0.5?
        elif vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        else:
            if clip:
                mask = np.ma.getmask(result)
                result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                     mask=mask)
            # ma division is very slow; we can take a shortcut
            resdat = result.data
            resdat -= vmin
            resdat /= (vmax - vmin)
            result = np.ma.array(resdat, mask=result.mask, copy=False)
        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        (vmin,), _ = self.process_value(self.vmin)
        (vmax,), _ = self.process_value(self.vmax)

        if cbook.iterable(value):
            val = np.ma.asarray(value)
            return vmin + val * (vmax - vmin)
        else:
            return vmin + value * (vmax - vmin)

    def autoscale(self, A):
        """
        Set *vmin*, *vmax* to min, max of *A*.
        """
        A = np.asanyarray(A)
        self.vmin = A.min()
        self.vmax = A.max()

    def autoscale_None(self, A):
        """autoscale only None-valued vmin or vmax."""
        A = np.asanyarray(A)
        if self.vmin is None and A.size:
            self.vmin = A.min()
        if self.vmax is None and A.size:
            self.vmax = A.max()

    def scaled(self):
        'return true if vmin and vmax set'
        return (self.vmin is not None and self.vmax is not None)

class PowerNorm(Normalize):
    """
    Normalize a given value to the ``[0, 1]`` interval with a power-law
    scaling. This will clip any negative data points to 0.
    """
    def __init__(self, gamma, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self, vmin, vmax, clip)
        self.gamma = gamma

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        gamma = self.gamma
        vmin, vmax = self.vmin, self.vmax
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            result.fill(0)
        else:
            res_mask = result.data < 0
            if clip:
                mask = np.ma.getmask(result)
                result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                     mask=mask)
            resdat = result.data
            resdat -= vmin
            np.power(resdat, gamma, resdat)
            resdat /= (vmax - vmin) ** gamma

            result = np.ma.array(resdat, mask=result.mask, copy=False)
            result[res_mask] = 0
        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        gamma = self.gamma
        vmin, vmax = self.vmin, self.vmax

        if cbook.iterable(value):
            val = np.ma.asarray(value)
            return np.ma.power(val, 1. / gamma) * (vmax - vmin) + vmin
        else:
            return pow(value, 1. / gamma) * (vmax - vmin) + vmin

    def autoscale(self, A):
        """
        Set *vmin*, *vmax* to min, max of *A*.
        """
        self.vmin = np.ma.min(A)
        if self.vmin < 0:
            self.vmin = 0
            warnings.warn("Power-law scaling on negative values is "
                          "ill-defined, clamping to 0.")
        self.vmax = np.ma.max(A)

    def autoscale_None(self, A):
        """autoscale only None-valued vmin or vmax."""
        A = np.asanyarray(A)
        if self.vmin is None and A.size:
            self.vmin = A.min()
            if self.vmin < 0:
                self.vmin = 0
                warnings.warn("Power-law scaling on negative values is "
                              "ill-defined, clamping to 0.")
        if self.vmax is None and A.size:
            self.vmax = A.max()