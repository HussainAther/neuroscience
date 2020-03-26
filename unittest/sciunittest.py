import math
import numpy as np
import sciunit

from datetime import datetime, timedelta
from math import pi, sqrt, sin, cos, tan, atan
from sciunit.capabilities import ProducesNumber

class ProducesOrbitalPosition(sciunit.Capability):
    """
    A model 'capability', i.e. a collection of methods that a test is allowed to invoke on a model.
    These methods are unimplemented by design, and the model must implement them.
    """
    
    def get_position(self, t: datetime) -> tuple:
        """
        Produce an orbital position from a time point
        in polar coordinates.
        
        Args:
            t (datetime): The time point to examine, relative to perihelion
        Returns:
            tuple: A pair of (r, theta) coordinates in the oribtal plane
        """
        raise NotImplementedError("")
        
    @property
    def perihelion(self) -> datetime:
        """
        Return the time of last perihelion.
        """
        raise NotImplementedError("")
    
    @property
    def period(self) -> float:
        """Return the period of the orbit"""
        raise NotImplementedError("")
    
    @property
    def eccentricity(self) -> float:
        """Return the eccentricity of the orbit"""
        raise NotImplementedError("")
    
    def get_x_y(self, t: datetime) -> tuple:
        """Produce an orbital position from a time point, but in cartesian coordinates.
        This method does not require a model-specific implementation.
        Thus, a generic implementation can be provided in advance."""
        r, theta = self.get_position(t)
        x, y = r*cos(theta), r*sin(theta)
        return x, y
