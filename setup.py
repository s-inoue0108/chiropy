from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="chiropy",
    version="1.0.0",
    description="CHIROpy is a Gaussian binding tool for analyze chiroptical properties.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Shota Inoue",
    author_email="inoue.shota@st.kitasato-u.ac.jp",
    entry_points={
        "console_scripts": [
            "chiropy = chiropy.main:main",
        ],
    },
    include_package_data=True,
    packages=find_packages(),
    install_requires=open("requirements.txt").read().splitlines(),
    python_requires=">=3.10",
)

