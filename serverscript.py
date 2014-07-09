from slab.instruments import InstrumentManager
from os import getcwd
print getcwd()
im=InstrumentManager(r'.\instruments.cfg', server=True)
