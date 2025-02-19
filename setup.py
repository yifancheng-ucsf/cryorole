from setuptools import setup, find_packages

setup(
    name="cryorole",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "scipy",
        "mpi4py"
    ],
    entry_points={
        'console_scripts': [
            'orientation_analysis=cryorole.OrientationAnalysis:main',
            'landscape_projection=cryorole.LandscapeProjection:main',
            'point_select=cryorole.PointSelect:main',
            'particle_backtrack=cryorole.ParticleBacktrack:main'
        ]
    },
    author="Chengmin Li",
    author_email="chengmin.li@ucsf.edu",
    description="A suite of tools to analyze and visualize orientation landscapes.",
    url="https://github.com/yifancheng-ucsf/cryorole",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License"
    ],
    python_requires='>=3.6',
)

