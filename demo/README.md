DEMOS
=====
These demonstrations are intended to highlight the pieces of the
`WholeBrain_MVPA` system and how they interoperate, so you will be able to
use the tools in the service of your own analysis. So far there are two demos.

demo00
------
This demo starts from nothing and completely generates data and metadata
objects for use with `WholeBrain_MVPA`. I recommend walking carefully
through this demo to get a clear sense for how to specify and run an
analysis. This demo does everything ``the hard way'' in an attempt to
foster intuitions about the system. However, certain parts of the
standard workflow have been abstracted into utility programs. This demo
generally does not use those.

demo01
------
This demo is the compliment to the first: it depends on you starting
with a particular dataset with well-formed metadata. At present, this
dataset is not publicly available and so this tutorial is not fully
functional for the public at large. The objective of this demo is a
walkthrough of a standard workflow, including operating the
`condortools` utilities for setting up jobs for use on a High Throughput
Computing distributed cluster. This demo builds on low level intuitions
that demo00 was written to foster.
