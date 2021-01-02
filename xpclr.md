I am going to use an approach called XPCLR to perform a genomewide scan for selection.

1. Download the code:
```
git clone https://github.com/hardingnj/xpclr.git
```
2. Load python module:
```
module load python
```
3. Create a virtual environment:
```
virtualenv /home/$USER/xpclr_venv
```
4. Activate your virtual environment:
```
source /home/$USER/xpclr_venv/bin/activate
```
5. Install the dependencies:

The package come with a file "requirements.txt" that contains the list of dependencies. Usually, the installation should work by running "pip install -r requirements.txt". However, in this case, numpy is required for installing "scikit-allel". Therefore, "numpy" should appear first on the list. You can just edit the file and change the order of the packages or install numpy first.

```
(xpclr_venv) [~]$ cd xpclr
(xpclr_venv) [xpclr]$ cat requirements.txt
scikit-allel>=1.2
numpy
pandas
scipy
h5py
zarr
```

To install the dependencies, I did run the following:

```
pip install numpy --no-index
pip install -r requirements.txt
```
6. Install the package "xpclr":
```
(xpclr_venv) [~]$ cd xpclr
(xpclr_venv) [~]$ python setup.py install
(xpclr_venv) [~]$ deactivate
```
7. Test the installation
```
[~]$ module load python
[~]$ source /home/$USER/xpclr_venv/bin/activate
(xpclr_venv) [~]$ xpclr --help
```
