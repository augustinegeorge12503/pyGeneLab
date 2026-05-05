-e git+ssh://git@github.com/augustinegeorge12503/pyGeneLab.git@main#egg=pyGeneLab

# first time
python -m pip install -r requirements.txt

# whenever github library updates
python -m pip install --upgrade --force-reinstall -r requirements.txt
python -m pip install --upgrade --no-deps --force-reinstall git+ssh://git@github.com/augustinegeorge12503/pyGeneLab.git
