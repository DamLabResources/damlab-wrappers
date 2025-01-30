PREFIX = $(PWD)/venv

DACCORD_URL = https://github.com/gt1/daccord/releases/download/0.0.10-release-20170526170720/daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu.tar.gz

install: $(PREFIX) $(PREFIX)/bin/daccord $(PREFIX)/bin/strainline.sh

$(PREFIX):
	mamba env create -p $(PREFIX) -f environment.yaml

$(PREFIX)/bin/daccord:
	mkdir -p $(PREFIX)/daccord
	wget -qO- $(DACCORD_URL) | tar xvz -C $(PREFIX)/
	ln -s $(PREFIX)/daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu/bin/daccord $(PREFIX)/bin/daccord

$(PREFIX)/bin/strainline.sh:
	git clone git@github.com:HaploKit/Strainline.git $(PREFIX)/strainline
	ln -s $(PREFIX)/strainline/src/strainline.sh $(PREFIX)/bin/strainline.sh
