# Test for fast_diffusion with a constant concentration

This is a test for the `fast_diffusion` in mode `FAST_DIFFUSION_RESERVOIR`.

There are three sets of simulations:

* To verify that works with normal singles.
* To verify that works with digits that do not compete.
* To verify that works with digits that compete for same binding sites.

Results are shown here:

![](plot0.svg)
![](plot1.svg)
![](plot2.svg)

The results can be reproduced by going to `test_reservoir` and running:

``` bash
bash run_sims.sh runs?
python plot.py

# If you want to delete all generated files:
bash cleanup_test.sh runs?
```


