from setuptools import setup, find_packages

setup(
    name='paint_atmospheres',
    version='0.1.0',
    packages=find_packages(),
    package_dir={'paint_atmospheres': 'pa'},
    entry_points={'console_scripts': [
    	'calc_limbdark=pa.calc_limbdark:run',
    	'calc_star=pa.calc_star:run',
    	'calc_spectra=pa.calc_spectra:run',
    	'plot_spectrum=pa.calc_spectrum:run',
        'plot_Imu=pa.opt.plot_Imu:run',
        'convert_spectrum=pa.opt.convert_spectrum:run'
    	]}
)
