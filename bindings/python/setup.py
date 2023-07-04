import distutils.core

distutils.core.setup(
	name = 'kww-python',
	author = 'Joachim Wuttke',
	author_email = 'j.wuttke@fz-juelich.de',
	url = 'http://apps.jcns.fz-juelich.de/kww',
	description = 'Computes the Kohlrausch-Williams-Watts (Fourier-Laplace transform of the stretched exponential function)',
	license = 'GPLv3',
	ext_modules = [
		distutils.core.Extension(
			'_kww',
			['kww.c', 'kww.i'],
			extra_compile_args = ['--std=c99']
		)
	],
	py_modules = [
		'kww'
	],
)
