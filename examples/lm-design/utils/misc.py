# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
import math
import sys
from contextlib import AbstractContextManager, contextmanager
import inspect
import logging
import os
from pathlib import Path
import pickle
import socket
import time

logger = logging.getLogger(__name__)



class TimerContext(AbstractContextManager):
    """
    Used to measure time and log if time threshold passed.
    Usage:
    with TimerContext(tag="disk usage", theshold=10.0) as tc:
       some_long_operation()  # will log warning in case this takes more than 10 seconds
    print(tc.elapsed) # elapsed data-member contains elapsed time in secs
    """
    def __init__(self, tag="N/A", threshold=None):
        super().__init__()
        self.tag = tag
        self.threshold = threshold

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.elapsed = (time.time() - self.start)
        self.elapsed_rounded = round(self.elapsed, 2)
        if self.threshold is not None and self.elapsed > self.threshold:
            logger.warning(f"Timer passed threshold: {self.tag} -- {self.elapsed_rounded}")
