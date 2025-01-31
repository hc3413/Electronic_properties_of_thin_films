import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime
import pathlib
from pathlib import Path
import math
import scipy
from scipy.optimize import fsolve
from scipy.stats import linregress
from pptx import Presentation
from pptx.util import Inches, Pt, Cm
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.oxml.ns import qn
from pptx.oxml import parse_xml
from dataclasses import dataclass, field
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import scipy.signal
from scipy.stats import zscore
import scienceplots


# Enable LaTeX for text rendering
plt.rcParams['text.usetex'] = True

plt.style.use(['science', 'grid'])








