import math
import numpy as np
import sciunit

from datetime import datetime, timedelta
from math import pi, sqrt, sin, cos, tan, atan
from sciunit.capabilities import ProducesNumber

"""
Testing SciUnit's abilities.
"""

class ProducesOrbitalPosition(sciunit.Capability):
    """
    A model "capability", i.e. a collection of methods that a test is allowed to invoke on a model.
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
        """
        Return the period of the orbit.
        """
        raise NotImplementedError("")
    
    @property
    def eccentricity(self) -> float:
        """
        Return the eccentricity of the orbit.
        """
        raise NotImplementedError("")
    
    def get_x_y(self, t: datetime) -> tuple:
        """
        Produce an orbital position from a time point, but in cartesian coordinates.
        This method does not require a model-specific implementation.
        Thus, a generic implementation can be provided in advance.
        """
        r, theta = self.get_position(t)
        x, y = r*cos(theta), r*sin(theta)
        return x, y

class PositionTest(sciunit.Test):
    """
    A test of a planetary position at some specified time.
    """
    # This test can only operate on models that implement
    # the "ProducesOrbitalPosition" capability.
    required_capabilities = (ProducesOrbitalPosition,)
    score_type = BooleanScore # This test"s "judge" method will return a BooleanScore.
    
    def generate_prediction(self, model):
        """
        Generate a prediction from a model.
        """
        t = self.observation["t"] # Get the time point from the test"s observation
        x, y = model.get_x_y(t) # Get the predicted x, y coordinates from the model
        return {"t": t, "x": x, "y": y} # Roll this into a model prediction dictionary
        
    def compute_score(self, observation, prediction):
        """
        Compute a test score based on the agreement between
        the observation (data) and prediction (model).
        """
        # Compare observation and prediction to get an error measure.
        delta_x = observation["x"] - prediction["x"]
        delta_y = observation["y"] - prediction["y"]
        error = np.sqrt(delta_x**2 + delta_y**2)
        passing = bool(error < 1e5*pq.kilometer) # Turn this into a True/False score.
        score = self.score_type(passing) # Create a sciunit.Score object.
        score.set_raw(error) # Add some information about how this score was obtained.
        score.description = ("Passing score if the prediction is "
                             "within < 100,000 km of the observation") # Describe the scoring logic.
        return score
