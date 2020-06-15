from setuptools import setup, find_packages

setup(
    name='paint_atmospheres',
    version='1.1.1',
    packages=find_packages(),
    package_dir={'paint_atmospheres': 'pa'},
    entry_points={'console_scripts': [
    	'calc_limbdark=pa.calc_limbdark:run',
    	'calc_star=pa.calc_star:run',
    	'calc_spectra=pa.calc_spectra:run'
    	]}
)
