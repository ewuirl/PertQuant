#!/usr/local/opt/python@3.7/bin/python3.7
# -*- coding: utf-8 -*-
import re
import sys
from PertQuant.analyzeFASTQ.NanoPlot_hist import hist_main
# Place this file in usr/local/bin
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(hist_main())