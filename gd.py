import numpy as np

class GradientDescent:
	def __init__(self, f, df, alpha = 0.01, eps = 1e-6):
		self.f = f
		self.df = df
		self.alpha = alpha #変化率:alpha
		self.eps = eps #誤差:eps
		self.path = None

	def solve(self, init):
		x = init #初期条件
		path = []
		grad = self.df(x) # fの勾配を返す	
		path.append(x)
		print(x)
		while (grad**2).sum() > self.eps**2: # ▽fが、ある一定まで小さくなったら停止(他の停止条件も入れた方がいいのでは？)
			x = x - self.alpha * grad
			grad = self.df(x)
			path.append(x)
		self.path__ = np.array(path)
		self.x__ = x
		self.opt__ = self.f(x)
