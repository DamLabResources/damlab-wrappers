from setuptools import setup, find_packages

setup(
    name='multiqc_damlab',
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
            'dorado = multiqc_damlab.modules.dorado:DoradoModule',
            'strainline = multiqc_damlab.modules.strainline:StrainlineModule',
        ],
        'multiqc.cli_options.v1': [
            'debug_option = multiqc_damlab.cli:debug_option'
        ],
    },
) 