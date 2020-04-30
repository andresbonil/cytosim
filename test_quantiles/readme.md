# Test the quantile function

Test that the quantile function does what we expect.

```bash
# To run the test
report fiber:quantiles hand=binder_hand>report_quantiles.txt
report single:attached>single_attached.txt
python do_test.py
```