# Forest Cover Map Segmentation

## Python development environment setup  
```shell
sudo apt-get update
sudo apt-get install python3-venv
sudo apt-get install python3-dev
sudo apt-get install -y python3-opencv
```


## GDAL setup
```shell
sudo add-apt-repository ppa:ubuntugis/ppa && sudo apt-get update
sudo apt-get update
sudo apt-get install gdal-bin
sudo apt-get install libgdal-dev
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal
```

## Python virtual environment
### Virtual environment creation
```shell
python -m venv fms-env
```
### Virtual environment activation
```shell
source ./fms-env/bin/activate
```

### Python package dependencies installation 
```shell
pip install --upgrade pip
pip install -r requirements.txt
```