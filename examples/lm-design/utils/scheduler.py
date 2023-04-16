# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from typing import Any, Optional, Union
from numbers import Number

import torch

SCHEDULER_REPOSITORY = {}

class Scheduler:
    def __init__(
        self, optimizer: Optional[Any], scheduler: str, initial: float, **schedulerspec
    ):
        if optimizer is None:
            # Dummy optimizer to wrap with lr_scheduler
            dummy = torch.tensor([], requires_grad=True)
            optimizer = torch.optim.SGD([dummy], lr=initial)
            self.dummy_optimizer = optimizer
        else:
            self.dummy_optimizer = None
        self.name = scheduler
        self._scheduler = getattr(torch.optim.lr_scheduler, self.name)(
            optimizer, **schedulerspec
        )
        self.initial = initial

    def step(self):
        if self.dummy_optimizer:
            self.dummy_optimizer.step()  # Avoid UserWarning about step order
        self._scheduler.step()

    def __call__(self):
        return self._scheduler.get_last_lr()[0]

SchedulerSpecDict = dict
SchedulerSpec = Union[float, str, SchedulerSpecDict]

class ConstantSchedule():
    """
    Schedule that returns a constant value, defined at init.
    Pickle-able, as opposed to `lambda: value`.
    """
    def __init__(self, value: Number):
        self.value = value
    def __call__(self):
        return self.value

def to_scheduler(sspec: SchedulerSpec):
    """ Helper for initializing schedulers from config. """
    if isinstance(sspec, str):
        sspec = SCHEDULER_REPOSITORY[sspec]
        optimizer = None
    elif isinstance(sspec, tuple):
        # There is an actual optimizer to wrap
        optimizer, sspec = sspec
    else:
        optimizer = None
    # Initialize real Scheduler, or dummy lambda
    if isinstance(sspec, Number):
        return ConstantSchedule(sspec)
    else:
        return Scheduler(optimizer, **sspec)

def set_scheduler_repo(repo: dict):
    global SCHEDULER_REPOSITORY
    SCHEDULER_REPOSITORY = repo
