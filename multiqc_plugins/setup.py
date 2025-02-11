from setuptools import setup, find_packages

setup(
    name='damlabqc',
    version='0.1.0',
    author='DAMlab',
    author_email='your.email@example.com',
    description='MultiQC plugin for DAMlab wrappers',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'multiqc>=1.14',
        'click>=7.0'
    ],
    entry_points={
        'multiqc.modules.v1': [
            'dorado = damlabqc.modules.dorado:DoradoModule',
            'strainline = damlabqc.modules.strainline:StrainlineModule',
        ],
        'multiqc.cli_options.v1': [
            'debug_option = damlabqc.cli:debug_option'
        ],
    },
) 