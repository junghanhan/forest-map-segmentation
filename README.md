# Forest Cover Map Segmentation

Please refer to `FMS_tool_user_guide.pdf` file for more detailed information.

## 1. Python development libraries installation  
```shell
sudo apt-get update
sudo apt-get install python3-venv
sudo apt-get install python3-dev
sudo apt-get install -y python3-opencv
```


## 2. GDAL development libraries installation (version 2.4.2)
```shell
sudo add-apt-repository ppa:ubuntugis/ppa && sudo apt-get update
sudo apt-get update
sudo apt-get install gdal-bin
sudo apt-get install libgdal-dev
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal
```

## 3. Cloning the FMS tool git repository into the host
```shell
git clone https://github.com/sigmabrains/forest-map-segmentation.git
```

## 4. Python virtual environment setup
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
pip install setuptools==57.5.0
pip install --upgrade pip
pip install -r requirements.txt
```

## 5. Running the FMS Tool
```shell
python src/main.py input/093E15e_1976_D_1.gtiff
```

You can process multiple GeoTiff files. They will be processed in parallel.

```shell
usage: main.py input_gtiff_file1 input_gtiff_file2 input_gtiff_file3 ...
```