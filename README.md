This is a demo for designing and implementing an LPV controller.

The conditions for synthesizing the LPV controller are from the following paper. <br>
Apkarian, P., & Adams, R. J. (1998). Advanced gain-scheduling techniques for uncertain systems. Control Systems Technology, IEEE Transactions on, 6(1), 21–32.

The numerical example used in the demo is from the following paper. <br>
Wu, F., Hua, X., Packard, A., & Becker, G. (1996). Induced L2 norm nontrol for LPV systems with bounded parameter variation rates, 998, 983–998.

LMILab is used to solve the LMI conditions but YALMIP plus other solvers such as Sedumi and SDPT3 is recommended due to much faster computation. 
