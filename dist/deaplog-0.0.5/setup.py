from setuptools import setup, find_packages

setup(
    name="deaplog",
    version="0.0.5",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "sympy",
        "fisher",
    ],
    python_requires=">=3.7",
) 