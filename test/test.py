import pyGBCuda as pylatt

quad = pylatt.Quad("qh1g2c30a", 0.268, -0.641957314648, NKick=100)
quad.run_qsympass4(100, 1000)