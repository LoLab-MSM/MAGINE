import pysb.macros as macros
from pysb import *


def catalyze(enz, sub, product, klist):
    return macros.catalyze(enz, 'bf', sub, 'bf', product, klist)


Model()
Monomer("hsa317", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa51807", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8793", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa3725", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8797", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8794", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8795", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa27113", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa2021", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa153090", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5783", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa7124", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1509", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa10039", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa10038", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa9131", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa4616", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8722", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa3845", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa142", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa143", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa112714", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa10912", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1513", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1512", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1515", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1514", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1519", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa84790", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa7132", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1147", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8717", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa355", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5894", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("cpdC00076", ["bf"])
Monomer("hsa9020", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa3551", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8772", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa572", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa356", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa4803", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5366", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa468", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa113457", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa27429", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1522", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa4000", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa4001", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("cpdC00319", ["bf"])
Monomer("hsa5551", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa2353", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa7186", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa7185", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8517", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa54205", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa4792", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa4790", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa208", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5970", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa207", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa100506742", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1647", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa23533", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1439", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1649", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1508", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8503", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa598", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa597", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa596", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa3562", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa3563", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa841", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa840", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa843", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa842", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1676", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1677", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1520", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa4893", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1075", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("cpdC05981", ["bf"])
Monomer("hsa472", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa63970", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5170", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa3002", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5595", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5594", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa2081", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5290", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5291", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5293", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5294", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa581", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5296", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa578", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa7278", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa4170", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa7277", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("cpdC00027", ["bf"])
Monomer("hsa7846", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa10018", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa10376", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa55367", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa3708", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa3709", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa10000", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa56616", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa84823", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa6709", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa6708", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa332", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa331", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa330", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa824", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa823", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8837", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1616", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5605", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5604", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5602", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5601", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8743", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa60", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa4914", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa637", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa9451", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa3710", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa79861", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5295", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa329", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5599", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa835", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa836", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa5414", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1521", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa839", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa7157", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa4217", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa1965", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8739", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa3265", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa8737", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Monomer("hsa71", ["gene", "state", "bf", ],
        {"state": ["I", "A"], "gene": ["on", "off", "protein"]})
Parameter("hsa317_0", 1)
Parameter("hsa51807_0", 1)
Parameter("hsa8793_0", 1)
Parameter("hsa3725_0", 1)
Parameter("hsa8797_0", 1)
Parameter("hsa8794_0", 1)
Parameter("hsa8795_0", 1)
Parameter("hsa27113_0", 1)
Parameter("hsa2021_0", 1)
Parameter("hsa153090_0", 1)
Parameter("hsa5783_0", 1)
Parameter("hsa7124_0", 1)
Parameter("hsa1509_0", 1)
Parameter("hsa10039_0", 1)
Parameter("hsa10038_0", 1)
Parameter("hsa9131_0", 1)
Parameter("hsa4616_0", 1)
Parameter("hsa8722_0", 1)
Parameter("hsa3845_0", 1)
Parameter("hsa142_0", 1)
Parameter("hsa143_0", 1)
Parameter("hsa112714_0", 1)
Parameter("hsa10912_0", 1)
Parameter("hsa1513_0", 1)
Parameter("hsa1512_0", 1)
Parameter("hsa1515_0", 1)
Parameter("hsa1514_0", 1)
Parameter("hsa1519_0", 1)
Parameter("hsa84790_0", 1)
Parameter("hsa7132_0", 1)
Parameter("hsa1147_0", 1)
Parameter("hsa8717_0", 1)
Parameter("hsa355_0", 1)
Parameter("hsa5894_0", 1)
Parameter("cpdC00076_0", 1)
Parameter("hsa9020_0", 1)
Parameter("hsa3551_0", 1)
Parameter("hsa8772_0", 1)
Parameter("hsa572_0", 1)
Parameter("hsa356_0", 1)
Parameter("hsa4803_0", 1)
Parameter("hsa5366_0", 1)
Parameter("hsa468_0", 1)
Parameter("hsa113457_0", 1)
Parameter("hsa27429_0", 1)
Parameter("hsa1522_0", 1)
Parameter("hsa4000_0", 1)
Parameter("hsa4001_0", 1)
Parameter("cpdC00319_0", 1)
Parameter("hsa5551_0", 1)
Parameter("hsa2353_0", 1)
Parameter("hsa7186_0", 1)
Parameter("hsa7185_0", 1)
Parameter("hsa8517_0", 1)
Parameter("hsa54205_0", 1)
Parameter("hsa4792_0", 1)
Parameter("hsa4790_0", 1)
Parameter("hsa208_0", 1)
Parameter("hsa5970_0", 1)
Parameter("hsa207_0", 1)
Parameter("hsa100506742_0", 1)
Parameter("hsa1647_0", 1)
Parameter("hsa23533_0", 1)
Parameter("hsa1439_0", 1)
Parameter("hsa1649_0", 1)
Parameter("hsa1508_0", 1)
Parameter("hsa8503_0", 1)
Parameter("hsa598_0", 1)
Parameter("hsa597_0", 1)
Parameter("hsa596_0", 1)
Parameter("hsa3562_0", 1)
Parameter("hsa3563_0", 1)
Parameter("hsa841_0", 1)
Parameter("hsa840_0", 1)
Parameter("hsa843_0", 1)
Parameter("hsa842_0", 1)
Parameter("hsa1676_0", 1)
Parameter("hsa1677_0", 1)
Parameter("hsa1520_0", 1)
Parameter("hsa4893_0", 1)
Parameter("hsa1075_0", 1)
Parameter("cpdC05981_0", 1)
Parameter("hsa472_0", 1)
Parameter("hsa63970_0", 1)
Parameter("hsa5170_0", 1)
Parameter("hsa3002_0", 1)
Parameter("hsa5595_0", 1)
Parameter("hsa5594_0", 1)
Parameter("hsa2081_0", 1)
Parameter("hsa5290_0", 1)
Parameter("hsa5291_0", 1)
Parameter("hsa5293_0", 1)
Parameter("hsa5294_0", 1)
Parameter("hsa581_0", 1)
Parameter("hsa5296_0", 1)
Parameter("hsa578_0", 1)
Parameter("hsa7278_0", 1)
Parameter("hsa4170_0", 1)
Parameter("hsa7277_0", 1)
Parameter("cpdC00027_0", 1)
Parameter("hsa7846_0", 1)
Parameter("hsa10018_0", 1)
Parameter("hsa10376_0", 1)
Parameter("hsa55367_0", 1)
Parameter("hsa3708_0", 1)
Parameter("hsa3709_0", 1)
Parameter("hsa10000_0", 1)
Parameter("hsa56616_0", 1)
Parameter("hsa84823_0", 1)
Parameter("hsa6709_0", 1)
Parameter("hsa6708_0", 1)
Parameter("hsa332_0", 1)
Parameter("hsa331_0", 1)
Parameter("hsa330_0", 1)
Parameter("hsa824_0", 1)
Parameter("hsa823_0", 1)
Parameter("hsa8837_0", 1)
Parameter("hsa1616_0", 1)
Parameter("hsa5605_0", 1)
Parameter("hsa5604_0", 1)
Parameter("hsa5602_0", 1)
Parameter("hsa5601_0", 1)
Parameter("hsa8743_0", 1)
Parameter("hsa60_0", 1)
Parameter("hsa4914_0", 1)
Parameter("hsa637_0", 1)
Parameter("hsa9451_0", 1)
Parameter("hsa3710_0", 1)
Parameter("hsa79861_0", 1)
Parameter("hsa5295_0", 1)
Parameter("hsa329_0", 1)
Parameter("hsa5599_0", 1)
Parameter("hsa835_0", 1)
Parameter("hsa836_0", 1)
Parameter("hsa5414_0", 1)
Parameter("hsa1521_0", 1)
Parameter("hsa839_0", 1)
Parameter("hsa7157_0", 1)
Parameter("hsa4217_0", 1)
Parameter("hsa1965_0", 1)
Parameter("hsa8739_0", 1)
Parameter("hsa3265_0", 1)
Parameter("hsa8737_0", 1)
Parameter("hsa71_0", 1)
Parameter("kf0", 1)
Parameter("kr0", 1)
Parameter("kc0", 1)
Parameter("kf1", 1)
Parameter("kr1", 1)
Parameter("kc1", 1)
Parameter("kf2", 1)
Parameter("kr2", 1)
Parameter("kc2", 1)
Parameter("k_synth2", 1)
Parameter("kf3", 1)
Parameter("kr3", 1)
Parameter("kc3", 1)
Parameter("k_synth3", 1)
Parameter("kf4", 1)
Parameter("kr4", 1)
Parameter("kc4", 1)
Parameter("k_synth4", 1)
Parameter("kf5", 1)
Parameter("kr5", 1)
Parameter("kc5", 1)
Parameter("k_synth5", 1)
Parameter("kf6", 1)
Parameter("kr6", 1)
Parameter("kc6", 1)
Parameter("k_synth6", 1)
Parameter("kf7", 1)
Parameter("kr7", 1)
Parameter("kc7", 1)
Parameter("kf8", 1)
Parameter("kr8", 1)
Parameter("kc8", 1)
Parameter("kf9", 1)
Parameter("kr9", 1)
Parameter("kc9", 1)
Parameter("kf10", 1)
Parameter("kr10", 1)
Parameter("kc10", 1)
Parameter("kf11", 1)
Parameter("kr11", 1)
Parameter("kc11", 1)
Parameter("kf12", 1)
Parameter("kr12", 1)
Parameter("kc12", 1)
Parameter("kf13", 1)
Parameter("kr13", 1)
Parameter("kc13", 1)
Parameter("kf14", 1)
Parameter("kr14", 1)
Parameter("kc14", 1)
Parameter("kf15", 1)
Parameter("kr15", 1)
Parameter("kc15", 1)
Parameter("kf16", 1)
Parameter("kr16", 1)
Parameter("kc16", 1)
Parameter("kf17", 1)
Parameter("kr17", 1)
Parameter("kc17", 1)
Parameter("kf18", 1)
Parameter("kr18", 1)
Parameter("kc18", 1)
Parameter("kf19", 1)
Parameter("kr19", 1)
Parameter("kc19", 1)
Parameter("kf20", 1)
Parameter("kr20", 1)
Parameter("kc20", 1)
Parameter("kf21", 1)
Parameter("kr21", 1)
Parameter("kc21", 1)
Parameter("kf22", 1)
Parameter("kr22", 1)
Parameter("kc22", 1)
Parameter("kf23", 1)
Parameter("kr23", 1)
Parameter("kc23", 1)
Parameter("kf24", 1)
Parameter("kr24", 1)
Parameter("kc24", 1)
Parameter("kf25", 1)
Parameter("kr25", 1)
Parameter("kc25", 1)
Parameter("kf26", 1)
Parameter("kr26", 1)
Parameter("kc26", 1)
Parameter("kf27", 1)
Parameter("kr27", 1)
Parameter("kc27", 1)
Parameter("kf28", 1)
Parameter("kr28", 1)
Parameter("kc28", 1)
Parameter("kf29", 1)
Parameter("kr29", 1)
Parameter("kc29", 1)
Parameter("kf30", 1)
Parameter("kr30", 1)
Parameter("kc30", 1)
Parameter("kf31", 1)
Parameter("kr31", 1)
Parameter("kc31", 1)
Parameter("kf32", 1)
Parameter("kr32", 1)
Parameter("kc32", 1)
Parameter("kf33", 1)
Parameter("kr33", 1)
Parameter("kc33", 1)
Parameter("kf34", 1)
Parameter("kr34", 1)
Parameter("kc34", 1)
Parameter("kf35", 1)
Parameter("kr35", 1)
Parameter("kc35", 1)
Parameter("kf36", 1)
Parameter("kr36", 1)
Parameter("kc36", 1)
Parameter("kf37", 1)
Parameter("kr37", 1)
Parameter("kc37", 1)
Parameter("kf38", 1)
Parameter("kr38", 1)
Parameter("kc38", 1)
Parameter("kf39", 1)
Parameter("kr39", 1)
Parameter("kc39", 1)
Parameter("kf40", 1)
Parameter("kr40", 1)
Parameter("kc40", 1)
Parameter("kf41", 1)
Parameter("kr41", 1)
Parameter("kc41", 1)
Parameter("kf42", 1)
Parameter("kr42", 1)
Parameter("kc42", 1)
Parameter("kf43", 1)
Parameter("kr43", 1)
Parameter("kc43", 1)
Parameter("kf44", 1)
Parameter("kr44", 1)
Parameter("kc44", 1)
Parameter("kf45", 1)
Parameter("kr45", 1)
Parameter("kc45", 1)
Parameter("kf46", 1)
Parameter("kr46", 1)
Parameter("kc46", 1)
Parameter("kf47", 1)
Parameter("kr47", 1)
Parameter("kc47", 1)
Parameter("kf48", 1)
Parameter("kr48", 1)
Parameter("kc48", 1)
Parameter("kf49", 1)
Parameter("kr49", 1)
Parameter("kc49", 1)
Parameter("kf50", 1)
Parameter("kr50", 1)
Parameter("kc50", 1)
Parameter("kf51", 1)
Parameter("kr51", 1)
Parameter("kc51", 1)
Parameter("kf52", 1)
Parameter("kr52", 1)
Parameter("kc52", 1)
Parameter("kf53", 1)
Parameter("kr53", 1)
Parameter("kc53", 1)
Parameter("kf54", 1)
Parameter("kr54", 1)
Parameter("kc54", 1)
Parameter("kf55", 1)
Parameter("kr55", 1)
Parameter("kc55", 1)
Parameter("kf56", 1)
Parameter("kr56", 1)
Parameter("kc56", 1)
Parameter("kf57", 1)
Parameter("kr57", 1)
Parameter("kc57", 1)
Parameter("kf58", 1)
Parameter("kr58", 1)
Parameter("kc58", 1)
Parameter("kf59", 1)
Parameter("kr59", 1)
Parameter("kc59", 1)
Parameter("kf60", 1)
Parameter("kr60", 1)
Parameter("kc60", 1)
Parameter("kf61", 1)
Parameter("kr61", 1)
Parameter("kc61", 1)
Parameter("kf62", 1)
Parameter("kr62", 1)
Parameter("kc62", 1)
Parameter("kf63", 1)
Parameter("kr63", 1)
Parameter("kc63", 1)
Parameter("kf64", 1)
Parameter("kr64", 1)
Parameter("kc64", 1)
Parameter("kf65", 1)
Parameter("kr65", 1)
Parameter("kc65", 1)
Parameter("kf66", 1)
Parameter("kr66", 1)
Parameter("kc66", 1)
Parameter("kf67", 1)
Parameter("kr67", 1)
Parameter("kc67", 1)
Parameter("kf68", 1)
Parameter("kr68", 1)
Parameter("kc68", 1)
Parameter("kf69", 1)
Parameter("kr69", 1)
Parameter("kc69", 1)
Parameter("kf70", 1)
Parameter("kr70", 1)
Parameter("kc70", 1)
Parameter("kf71", 1)
Parameter("kr71", 1)
Parameter("kc71", 1)
Parameter("kf72", 1)
Parameter("kr72", 1)
Parameter("kc72", 1)
Parameter("kf73", 1)
Parameter("kr73", 1)
Parameter("kc73", 1)
Parameter("kf74", 1)
Parameter("kr74", 1)
Parameter("kc74", 1)
Parameter("kf75", 1)
Parameter("kr75", 1)
Parameter("kc75", 1)
Parameter("kf76", 1)
Parameter("kr76", 1)
Parameter("kc76", 1)
Parameter("kf77", 1)
Parameter("kr77", 1)
Parameter("kc77", 1)
Parameter("kf78", 1)
Parameter("kr78", 1)
Parameter("kc78", 1)
Parameter("kf79", 1)
Parameter("kr79", 1)
Parameter("kc79", 1)
Parameter("kf80", 1)
Parameter("kr80", 1)
Parameter("kc80", 1)
Parameter("kf81", 1)
Parameter("kr81", 1)
Parameter("kc81", 1)
Parameter("kf82", 1)
Parameter("kr82", 1)
Parameter("kc82", 1)
Parameter("k_synth82", 1)
Parameter("kf83", 1)
Parameter("kr83", 1)
Parameter("kc83", 1)
Parameter("kf84", 1)
Parameter("kr84", 1)
Parameter("kc84", 1)
Parameter("kf85", 1)
Parameter("kr85", 1)
Parameter("kc85", 1)
Parameter("kf86", 1)
Parameter("kr86", 1)
Parameter("kc86", 1)
Parameter("kf87", 1)
Parameter("kr87", 1)
Parameter("kc87", 1)
Parameter("kf88", 1)
Parameter("kr88", 1)
Parameter("kc88", 1)
Parameter("kf89", 1)
Parameter("kr89", 1)
Parameter("kc89", 1)
Parameter("kf90", 1)
Parameter("kr90", 1)
Parameter("kc90", 1)
Parameter("kf91", 1)
Parameter("kr91", 1)
Parameter("kc91", 1)
Parameter("kf92", 1)
Parameter("kr92", 1)
Parameter("kc92", 1)
Parameter("kf93", 1)
Parameter("kr93", 1)
Parameter("kc93", 1)
Parameter("kf94", 1)
Parameter("kr94", 1)
Parameter("kc94", 1)
Parameter("k_synth94", 1)
Parameter("kf95", 1)
Parameter("kr95", 1)
Parameter("kc95", 1)
Parameter("k_synth95", 1)
Parameter("kf96", 1)
Parameter("kr96", 1)
Parameter("kc96", 1)
Parameter("k_synth96", 1)
Parameter("kf97", 1)
Parameter("kr97", 1)
Parameter("kc97", 1)
Parameter("k_synth97", 1)
Parameter("kf98", 1)
Parameter("kr98", 1)
Parameter("kc98", 1)
Parameter("k_synth98", 1)
Parameter("kf99", 1)
Parameter("kr99", 1)
Parameter("kc99", 1)
Parameter("kf100", 1)
Parameter("kr100", 1)
Parameter("kc100", 1)
Parameter("kf101", 1)
Parameter("kr101", 1)
Parameter("kc101", 1)
Parameter("kf102", 1)
Parameter("kr102", 1)
Parameter("kf103", 1)
Parameter("kr103", 1)
Parameter("kc103", 1)
Parameter("kf106", 1)
Parameter("kr106", 1)
Parameter("kc106", 1)
Parameter("k_synth106", 1)
Parameter("kf107", 1)
Parameter("kr107", 1)
Parameter("kc107", 1)
Parameter("k_synth107", 1)
Parameter("kf108", 1)
Parameter("kr108", 1)
Parameter("kc108", 1)
Parameter("k_synth108", 1)
Parameter("kf109", 1)
Parameter("kr109", 1)
Parameter("kc109", 1)
Parameter("k_synth109", 1)
Parameter("kf110", 1)
Parameter("kr110", 1)
Parameter("kc110", 1)
Parameter("k_synth110", 1)
Parameter("kf111", 1)
Parameter("kr111", 1)
Parameter("kc111", 1)
Parameter("k_synth111", 1)
Parameter("kf112", 1)
Parameter("kr112", 1)
Parameter("kc112", 1)
Parameter("k_synth112", 1)
Parameter("kf113", 1)
Parameter("kr113", 1)
Parameter("kc113", 1)
Parameter("k_synth113", 1)
Parameter("kf114", 1)
Parameter("kr114", 1)
Parameter("kc114", 1)
Parameter("k_synth114", 1)
Parameter("kf115", 1)
Parameter("kr115", 1)
Parameter("kc115", 1)
Parameter("k_synth115", 1)
Parameter("kf116", 1)
Parameter("kr116", 1)
Parameter("kc116", 1)
Parameter("k_synth116", 1)
Parameter("kf117", 1)
Parameter("kr117", 1)
Parameter("kc117", 1)
Parameter("k_synth117", 1)
Parameter("kf118", 1)
Parameter("kr118", 1)
Parameter("kc118", 1)
Parameter("k_synth118", 1)
Parameter("kf119", 1)
Parameter("kr119", 1)
Parameter("kc119", 1)
Parameter("kf120", 1)
Parameter("kr120", 1)
Parameter("kc120", 1)
Parameter("kf121", 1)
Parameter("kr121", 1)
Parameter("kc121", 1)
Parameter("kf122", 1)
Parameter("kr122", 1)
Parameter("kc122", 1)
Parameter("kf123", 1)
Parameter("kr123", 1)
Parameter("kc123", 1)
Parameter("k_synth123", 1)
Parameter("kf124", 1)
Parameter("kr124", 1)
Parameter("kc124", 1)
Parameter("k_synth124", 1)
Parameter("kf125", 1)
Parameter("kr125", 1)
Parameter("kc125", 1)
Parameter("k_synth125", 1)
Parameter("kf126", 1)
Parameter("kr126", 1)
Parameter("kc126", 1)
Parameter("k_synth126", 1)
Parameter("kf127", 1)
Parameter("kr127", 1)
Parameter("kc127", 1)
Parameter("k_synth127", 1)
Parameter("kf128", 1)
Parameter("kr128", 1)
Parameter("kc128", 1)
Parameter("k_synth128", 1)
Parameter("kf129", 1)
Parameter("kr129", 1)
Parameter("kc129", 1)
Parameter("k_synth129", 1)
Parameter("kf130", 1)
Parameter("kr130", 1)
Parameter("kc130", 1)
Parameter("k_synth130", 1)
Parameter("kf131", 1)
Parameter("kr131", 1)
Parameter("kc131", 1)
Parameter("k_synth131", 1)
Parameter("kf132", 1)
Parameter("kr132", 1)
Parameter("kc132", 1)
Parameter("k_synth132", 1)
Parameter("kf133", 1)
Parameter("kr133", 1)
Parameter("kc133", 1)
Parameter("k_synth133", 1)
Parameter("kf134", 1)
Parameter("kr134", 1)
Parameter("kc134", 1)
Parameter("k_synth134", 1)
Parameter("kf135", 1)
Parameter("kr135", 1)
Parameter("kc135", 1)
Parameter("k_synth135", 1)
Parameter("kf136", 1)
Parameter("kr136", 1)
Parameter("kc136", 1)
Parameter("kf137", 1)
Parameter("kr137", 1)
Parameter("kc137", 1)
Parameter("kf138", 1)
Parameter("kr138", 1)
Parameter("kc138", 1)
Parameter("kf139", 1)
Parameter("kr139", 1)
Parameter("kc139", 1)
Parameter("kf140", 1)
Parameter("kr140", 1)
Parameter("kc140", 1)
Parameter("kf141", 1)
Parameter("kr141", 1)
Parameter("kc141", 1)
Parameter("kf142", 1)
Parameter("kr142", 1)
Parameter("kc142", 1)
Parameter("kf143", 1)
Parameter("kr143", 1)
Parameter("kc143", 1)
Parameter("kf144", 1)
Parameter("kr144", 1)
Parameter("kc144", 1)
Parameter("kf145", 1)
Parameter("kr145", 1)
Parameter("kc145", 1)
Parameter("kf146", 1)
Parameter("kr146", 1)
Parameter("kc146", 1)
Parameter("kf147", 1)
Parameter("kr147", 1)
Parameter("kc147", 1)
Parameter("kf148", 1)
Parameter("kr148", 1)
Parameter("kc148", 1)
Parameter("kf149", 1)
Parameter("kr149", 1)
Parameter("kc149", 1)
Parameter("kf150", 1)
Parameter("kr150", 1)
Parameter("kc150", 1)
Parameter("kf151", 1)
Parameter("kr151", 1)
Parameter("kc151", 1)
Parameter("kf152", 1)
Parameter("kr152", 1)
Parameter("kc152", 1)
Parameter("kf153", 1)
Parameter("kr153", 1)
Parameter("kc153", 1)
Parameter("kf154", 1)
Parameter("kr154", 1)
Parameter("kc154", 1)
Parameter("kf155", 1)
Parameter("kr155", 1)
Parameter("kc155", 1)
Parameter("kf156", 1)
Parameter("kr156", 1)
Parameter("kc156", 1)
Parameter("kf157", 1)
Parameter("kr157", 1)
Parameter("kc157", 1)
Parameter("kf158", 1)
Parameter("kr158", 1)
Parameter("kc158", 1)
Parameter("kf159", 1)
Parameter("kr159", 1)
Parameter("kc159", 1)
Parameter("kf160", 1)
Parameter("kr160", 1)
Parameter("kc160", 1)
Parameter("kf161", 1)
Parameter("kr161", 1)
Parameter("kc161", 1)
Parameter("kf162", 1)
Parameter("kr162", 1)
Parameter("kc162", 1)
Parameter("kf163", 1)
Parameter("kr163", 1)
Parameter("kc163", 1)
Parameter("kf164", 1)
Parameter("kr164", 1)
Parameter("kc164", 1)
Parameter("kf165", 1)
Parameter("kr165", 1)
Parameter("kc165", 1)
Parameter("kf166", 1)
Parameter("kr166", 1)
Parameter("kc166", 1)
Parameter("kf167", 1)
Parameter("kr167", 1)
Parameter("kc167", 1)
Parameter("kf168", 1)
Parameter("kr168", 1)
Parameter("kc168", 1)
Parameter("kf169", 1)
Parameter("kr169", 1)
Parameter("kc169", 1)
Parameter("kf170", 1)
Parameter("kr170", 1)
Parameter("kc170", 1)
Parameter("kf171", 1)
Parameter("kr171", 1)
Parameter("kc171", 1)
Parameter("kf172", 1)
Parameter("kr172", 1)
Parameter("kc172", 1)
Parameter("kf173", 1)
Parameter("kr173", 1)
Parameter("kc173", 1)
Parameter("kf174", 1)
Parameter("kr174", 1)
Parameter("kc174", 1)
Parameter("kf175", 1)
Parameter("kr175", 1)
Parameter("kc175", 1)
Parameter("kf176", 1)
Parameter("kr176", 1)
Parameter("kc176", 1)
Parameter("kf177", 1)
Parameter("kr177", 1)
Parameter("kc177", 1)
Parameter("kf178", 1)
Parameter("kr178", 1)
Parameter("kc178", 1)
Parameter("kf179", 1)
Parameter("kr179", 1)
Parameter("kc179", 1)
Parameter("kf180", 1)
Parameter("kr180", 1)
Parameter("kc180", 1)
Parameter("kf181", 1)
Parameter("kr181", 1)
Parameter("kc181", 1)
Parameter("kf182", 1)
Parameter("kr182", 1)
Parameter("kc182", 1)
Parameter("kf183", 1)
Parameter("kr183", 1)
Parameter("kc183", 1)
Parameter("kf184", 1)
Parameter("kr184", 1)
Parameter("kc184", 1)
Parameter("kf185", 1)
Parameter("kr185", 1)
Parameter("kc185", 1)
Parameter("kf186", 1)
Parameter("kr186", 1)
Parameter("kc186", 1)
Parameter("kf187", 1)
Parameter("kr187", 1)
Parameter("kc187", 1)
Parameter("kf188", 1)
Parameter("kr188", 1)
Parameter("kc188", 1)
Parameter("kf189", 1)
Parameter("kr189", 1)
Parameter("kc189", 1)
Parameter("kf190", 1)
Parameter("kr190", 1)
Parameter("kc190", 1)
Parameter("kf191", 1)
Parameter("kr191", 1)
Parameter("kc191", 1)
Parameter("kf192", 1)
Parameter("kr192", 1)
Parameter("kc192", 1)
Parameter("kf193", 1)
Parameter("kr193", 1)
Parameter("kc193", 1)
Parameter("kf194", 1)
Parameter("kr194", 1)
Parameter("kc194", 1)
Parameter("kf195", 1)
Parameter("kr195", 1)
Parameter("kc195", 1)
Parameter("kf196", 1)
Parameter("kr196", 1)
Parameter("kc196", 1)
Parameter("kf197", 1)
Parameter("kr197", 1)
Parameter("kc197", 1)
Parameter("kf198", 1)
Parameter("kr198", 1)
Parameter("kc198", 1)
Parameter("kf199", 1)
Parameter("kr199", 1)
Parameter("kc199", 1)
Parameter("kf200", 1)
Parameter("kr200", 1)
Parameter("kc200", 1)
Parameter("kf201", 1)
Parameter("kr201", 1)
Parameter("kc201", 1)
Parameter("kf202", 1)
Parameter("kr202", 1)
Parameter("kc202", 1)
Parameter("kf203", 1)
Parameter("kr203", 1)
Parameter("kc203", 1)
Parameter("kf204", 1)
Parameter("kr204", 1)
Parameter("kc204", 1)
Parameter("kf205", 1)
Parameter("kr205", 1)
Parameter("kc205", 1)
Parameter("kf206", 1)
Parameter("kr206", 1)
Parameter("kc206", 1)
Parameter("kf207", 1)
Parameter("kr207", 1)
Parameter("kc207", 1)
Parameter("kf208", 1)
Parameter("kr208", 1)
Parameter("kc208", 1)
Parameter("kf209", 1)
Parameter("kr209", 1)
Parameter("kc209", 1)
Parameter("kf210", 1)
Parameter("kr210", 1)
Parameter("kc210", 1)
Parameter("kf211", 1)
Parameter("kr211", 1)
Parameter("kc211", 1)
Parameter("kf212", 1)
Parameter("kr212", 1)
Parameter("kc212", 1)
Parameter("kf213", 1)
Parameter("kr213", 1)
Parameter("kc213", 1)
Parameter("kf214", 1)
Parameter("kr214", 1)
Parameter("kc214", 1)
Parameter("kf215", 1)
Parameter("kr215", 1)
Parameter("kc215", 1)
Parameter("kf216", 1)
Parameter("kr216", 1)
Parameter("kc216", 1)
Parameter("kf217", 1)
Parameter("kr217", 1)
Parameter("kc217", 1)
Parameter("kf218", 1)
Parameter("kr218", 1)
Parameter("kc218", 1)
Parameter("kf219", 1)
Parameter("kr219", 1)
Parameter("kc219", 1)
Parameter("kf220", 1)
Parameter("kr220", 1)
Parameter("kc220", 1)
Parameter("kf221", 1)
Parameter("kr221", 1)
Parameter("kc221", 1)
Parameter("kf222", 1)
Parameter("kr222", 1)
Parameter("kc222", 1)
Parameter("kf223", 1)
Parameter("kr223", 1)
Parameter("kc223", 1)
Parameter("kf224", 1)
Parameter("kr224", 1)
Parameter("kc224", 1)
Parameter("kf225", 1)
Parameter("kr225", 1)
Parameter("kc225", 1)
Parameter("kf226", 1)
Parameter("kr226", 1)
Parameter("kc226", 1)
Parameter("kf227", 1)
Parameter("kr227", 1)
Parameter("kc227", 1)
Parameter("kf228", 1)
Parameter("kr228", 1)
Parameter("kc228", 1)
Parameter("kf229", 1)
Parameter("kr229", 1)
Parameter("kc229", 1)
Parameter("kf230", 1)
Parameter("kr230", 1)
Parameter("kc230", 1)
Parameter("kf231", 1)
Parameter("kr231", 1)
Parameter("kc231", 1)
Parameter("kf232", 1)
Parameter("kr232", 1)
Parameter("kc232", 1)
Parameter("kf233", 1)
Parameter("kr233", 1)
Parameter("kc233", 1)
Parameter("kf234", 1)
Parameter("kr234", 1)
Parameter("kc234", 1)
Parameter("kf235", 1)
Parameter("kr235", 1)
Parameter("kc235", 1)
Parameter("kf236", 1)
Parameter("kr236", 1)
Parameter("kc236", 1)
Parameter("kf237", 1)
Parameter("kr237", 1)
Parameter("kc237", 1)
Parameter("kf238", 1)
Parameter("kr238", 1)
Parameter("kc238", 1)
Parameter("kf239", 1)
Parameter("kr239", 1)
Parameter("kc239", 1)
Parameter("k_synth239", 1)
Parameter("kf240", 1)
Parameter("kr240", 1)
Parameter("kc240", 1)
Parameter("k_synth240", 1)
Parameter("kf241", 1)
Parameter("kr241", 1)
Parameter("kc241", 1)
Parameter("kf242", 1)
Parameter("kr242", 1)
Parameter("kc242", 1)
Parameter("kf243", 1)
Parameter("kr243", 1)
Parameter("kc243", 1)
Parameter("kf244", 1)
Parameter("kr244", 1)
Parameter("kc244", 1)
Parameter("kf245", 1)
Parameter("kr245", 1)
Parameter("kc245", 1)
Parameter("kf246", 1)
Parameter("kr246", 1)
Parameter("kc246", 1)
Parameter("kf247", 1)
Parameter("kr247", 1)
Parameter("kc247", 1)
Parameter("kf248", 1)
Parameter("kr248", 1)
Parameter("kc248", 1)
Parameter("kf249", 1)
Parameter("kr249", 1)
Parameter("kc249", 1)
Parameter("kf250", 1)
Parameter("kr250", 1)
Parameter("kc250", 1)
Parameter("kf251", 1)
Parameter("kr251", 1)
Parameter("kc251", 1)
Parameter("kf252", 1)
Parameter("kr252", 1)
Parameter("kc252", 1)
Parameter("kf253", 1)
Parameter("kr253", 1)
Parameter("kc253", 1)
Parameter("kf254", 1)
Parameter("kr254", 1)
Parameter("kc254", 1)
Parameter("kf255", 1)
Parameter("kr255", 1)
Parameter("kc255", 1)
Parameter("kf256", 1)
Parameter("kr256", 1)
Parameter("kc256", 1)
Parameter("kf257", 1)
Parameter("kr257", 1)
Parameter("kc257", 1)
Parameter("kf258", 1)
Parameter("kr258", 1)
Parameter("kc258", 1)
Parameter("kf259", 1)
Parameter("kr259", 1)
Parameter("kc259", 1)
Parameter("kf260", 1)
Parameter("kr260", 1)
Parameter("kc260", 1)
Parameter("kf261", 1)
Parameter("kr261", 1)
Parameter("kc261", 1)
Parameter("kf262", 1)
Parameter("kr262", 1)
Parameter("kc262", 1)
Parameter("kf263", 1)
Parameter("kr263", 1)
Parameter("kc263", 1)
Parameter("kf264", 1)
Parameter("kr264", 1)
Parameter("kc264", 1)
Parameter("kf265", 1)
Parameter("kr265", 1)
Parameter("kc265", 1)
Parameter("kf266", 1)
Parameter("kr266", 1)
Parameter("kc266", 1)
Parameter("kf267", 1)
Parameter("kr267", 1)
Parameter("kc267", 1)
Parameter("kf268", 1)
Parameter("kr268", 1)
Parameter("kc268", 1)
Parameter("kf269", 1)
Parameter("kr269", 1)
Parameter("kc269", 1)
Parameter("kf270", 1)
Parameter("kr270", 1)
Parameter("kc270", 1)
Parameter("kf271", 1)
Parameter("kr271", 1)
Parameter("kc271", 1)
Parameter("kf272", 1)
Parameter("kr272", 1)
Parameter("kc272", 1)
Parameter("kf273", 1)
Parameter("kr273", 1)
Parameter("kc273", 1)
Parameter("kf274", 1)
Parameter("kr274", 1)
Parameter("kc274", 1)
Parameter("kf275", 1)
Parameter("kr275", 1)
Parameter("kc275", 1)
Parameter("kf276", 1)
Parameter("kr276", 1)
Parameter("kc276", 1)
Parameter("kf277", 1)
Parameter("kr277", 1)
Parameter("kc277", 1)
Parameter("kf278", 1)
Parameter("kr278", 1)
Parameter("kc278", 1)
Parameter("kf279", 1)
Parameter("kr279", 1)
Parameter("kc279", 1)
Parameter("kf280", 1)
Parameter("kr280", 1)
Parameter("kc280", 1)
Parameter("kf281", 1)
Parameter("kr281", 1)
Parameter("kc281", 1)
Parameter("kf282", 1)
Parameter("kr282", 1)
Parameter("kc282", 1)
Parameter("kf283", 1)
Parameter("kr283", 1)
Parameter("kc283", 1)
Parameter("kf284", 1)
Parameter("kr284", 1)
Parameter("kc284", 1)
Parameter("kf285", 1)
Parameter("kr285", 1)
Parameter("kc285", 1)
Parameter("kf286", 1)
Parameter("kr286", 1)
Parameter("kc286", 1)
Parameter("kf287", 1)
Parameter("kr287", 1)
Parameter("kc287", 1)
Parameter("kf288", 1)
Parameter("kr288", 1)
Parameter("kc288", 1)
Parameter("kf289", 1)
Parameter("kr289", 1)
Parameter("kc289", 1)
Parameter("kf290", 1)
Parameter("kr290", 1)
Parameter("kc290", 1)
Parameter("kf291", 1)
Parameter("kr291", 1)
Parameter("kc291", 1)
Parameter("kf292", 1)
Parameter("kr292", 1)
Parameter("kc292", 1)
Parameter("kf293", 1)
Parameter("kr293", 1)
Parameter("kc293", 1)
Parameter("kf294", 1)
Parameter("kr294", 1)
Parameter("kc294", 1)
Parameter("kf295", 1)
Parameter("kr295", 1)
Parameter("kc295", 1)
Parameter("kf296", 1)
Parameter("kr296", 1)
Parameter("kc296", 1)
Parameter("kf297", 1)
Parameter("kr297", 1)
Parameter("kc297", 1)
Parameter("kf298", 1)
Parameter("kr298", 1)
Parameter("kc298", 1)
Parameter("kf299", 1)
Parameter("kr299", 1)
Parameter("kc299", 1)
Parameter("kf300", 1)
Parameter("kr300", 1)
Parameter("kc300", 1)
Parameter("kf301", 1)
Parameter("kr301", 1)
Parameter("kc301", 1)
Parameter("kf302", 1)
Parameter("kr302", 1)
Parameter("kc302", 1)
Parameter("kf303", 1)
Parameter("kr303", 1)
Parameter("kc303", 1)
Parameter("kf304", 1)
Parameter("kr304", 1)
Parameter("kc304", 1)
Parameter("kf305", 1)
Parameter("kr305", 1)
Parameter("kc305", 1)
Parameter("kf306", 1)
Parameter("kr306", 1)
Parameter("kc306", 1)
Parameter("kf307", 1)
Parameter("kr307", 1)
Parameter("kc307", 1)
Parameter("kf308", 1)
Parameter("kr308", 1)
Parameter("kc308", 1)
Parameter("kf309", 1)
Parameter("kr309", 1)
Parameter("kc309", 1)
Parameter("kf310", 1)
Parameter("kr310", 1)
Parameter("kc310", 1)
Parameter("kf311", 1)
Parameter("kr311", 1)
Parameter("kc311", 1)
Parameter("kf312", 1)
Parameter("kr312", 1)
Parameter("kc312", 1)
Parameter("kf313", 1)
Parameter("kr313", 1)
Parameter("kc313", 1)
Parameter("kf314", 1)
Parameter("kr314", 1)
Parameter("kc314", 1)
Parameter("kf315", 1)
Parameter("kr315", 1)
Parameter("kc315", 1)
Parameter("kf316", 1)
Parameter("kr316", 1)
Parameter("kc316", 1)
Parameter("kf317", 1)
Parameter("kr317", 1)
Parameter("kc317", 1)
Parameter("kf318", 1)
Parameter("kr318", 1)
Parameter("kc318", 1)
Parameter("kf319", 1)
Parameter("kr319", 1)
Parameter("kc319", 1)
Parameter("kf320", 1)
Parameter("kr320", 1)
Parameter("kc320", 1)
Parameter("kf321", 1)
Parameter("kr321", 1)
Parameter("kc321", 1)
Parameter("kf322", 1)
Parameter("kr322", 1)
Parameter("kc322", 1)
Parameter("kf323", 1)
Parameter("kr323", 1)
Parameter("kc323", 1)
Parameter("kf324", 1)
Parameter("kr324", 1)
Parameter("kc324", 1)
Parameter("kf325", 1)
Parameter("kr325", 1)
Parameter("kc325", 1)
Parameter("kf326", 1)
Parameter("kr326", 1)
Parameter("kc326", 1)
Parameter("kf327", 1)
Parameter("kr327", 1)
Parameter("kc327", 1)
Parameter("kf328", 1)
Parameter("kr328", 1)
Parameter("kc328", 1)
Parameter("kf329", 1)
Parameter("kr329", 1)
Parameter("kc329", 1)
Parameter("kf330", 1)
Parameter("kr330", 1)
Parameter("kc330", 1)
Parameter("kf331", 1)
Parameter("kr331", 1)
Parameter("kc331", 1)
Parameter("kf332", 1)
Parameter("kr332", 1)
Parameter("kc332", 1)
Parameter("kf333", 1)
Parameter("kr333", 1)
Parameter("kc333", 1)
Parameter("kf334", 1)
Parameter("kr334", 1)
Parameter("kc334", 1)
Parameter("kf335", 1)
Parameter("kr335", 1)
Parameter("kc335", 1)
Parameter("kf336", 1)
Parameter("kr336", 1)
Parameter("kc336", 1)
Parameter("kf337", 1)
Parameter("kr337", 1)
Parameter("kc337", 1)
Parameter("kf338", 1)
Parameter("kr338", 1)
Parameter("kc338", 1)
Parameter("kf339", 1)
Parameter("kr339", 1)
Parameter("kc339", 1)
Parameter("k_synth339", 1)
Parameter("kf340", 1)
Parameter("kr340", 1)
Parameter("kc340", 1)
Parameter("k_synth340", 1)
Parameter("kf341", 1)
Parameter("kr341", 1)
Parameter("kc341", 1)
Parameter("k_synth341", 1)
Parameter("kf342", 1)
Parameter("kr342", 1)
Parameter("kc342", 1)
Parameter("k_synth342", 1)
Parameter("kf343", 1)
Parameter("kr343", 1)
Parameter("kc343", 1)
Parameter("k_synth343", 1)
Parameter("kf344", 1)
Parameter("kr344", 1)
Parameter("kc344", 1)
Parameter("k_synth344", 1)
Parameter("kf345", 1)
Parameter("kr345", 1)
Parameter("kc345", 1)
Parameter("k_synth345", 1)
Parameter("kf346", 1)
Parameter("kr346", 1)
Parameter("kc346", 1)
Parameter("k_synth346", 1)
Parameter("kf347", 1)
Parameter("kr347", 1)
Parameter("kc347", 1)
Parameter("k_synth347", 1)
Parameter("kf348", 1)
Parameter("kr348", 1)
Parameter("kc348", 1)
Parameter("k_synth348", 1)
Parameter("kf349", 1)
Parameter("kr349", 1)
Parameter("kc349", 1)
Parameter("k_synth349", 1)
Parameter("kf350", 1)
Parameter("kr350", 1)
Parameter("kc350", 1)
Parameter("k_synth350", 1)
Parameter("kf351", 1)
Parameter("kr351", 1)
Parameter("kc351", 1)
Parameter("k_synth351", 1)
Parameter("kf352", 1)
Parameter("kr352", 1)
Parameter("kc352", 1)
Parameter("kf353", 1)
Parameter("kr353", 1)
Parameter("kc353", 1)
Parameter("kf354", 1)
Parameter("kr354", 1)
Parameter("kc354", 1)
Parameter("kf355", 1)
Parameter("kr355", 1)
Parameter("kc355", 1)
Parameter("kf356", 1)
Parameter("kr356", 1)
Parameter("kc356", 1)
Initial(hsa317(bf=None, gene="on", state="A"), hsa317_0)
Initial(hsa317(bf=None, gene="protein", state="A"), hsa317_0)
Initial(hsa317(bf=None, gene="on", state="I"), hsa317_0)
Initial(hsa317(bf=None, gene="off", state="I"), hsa317_0)
Initial(hsa51807(bf=None, gene="on", state="A"), hsa51807_0)
Initial(hsa51807(bf=None, gene="protein", state="A"), hsa51807_0)
Initial(hsa51807(bf=None, gene="on", state="I"), hsa51807_0)
Initial(hsa51807(bf=None, gene="off", state="I"), hsa51807_0)
Initial(hsa8793(bf=None, gene="on", state="A"), hsa8793_0)
Initial(hsa8793(bf=None, gene="protein", state="A"), hsa8793_0)
Initial(hsa8793(bf=None, gene="on", state="I"), hsa8793_0)
Initial(hsa8793(bf=None, gene="off", state="I"), hsa8793_0)
Initial(hsa3725(bf=None, gene="on", state="A"), hsa3725_0)
Initial(hsa3725(bf=None, gene="protein", state="A"), hsa3725_0)
Initial(hsa3725(bf=None, gene="on", state="I"), hsa3725_0)
Initial(hsa3725(bf=None, gene="off", state="I"), hsa3725_0)
Initial(hsa8797(bf=None, gene="on", state="A"), hsa8797_0)
Initial(hsa8797(bf=None, gene="protein", state="A"), hsa8797_0)
Initial(hsa8797(bf=None, gene="on", state="I"), hsa8797_0)
Initial(hsa8797(bf=None, gene="off", state="I"), hsa8797_0)
Initial(hsa8794(bf=None, gene="on", state="A"), hsa8794_0)
Initial(hsa8794(bf=None, gene="protein", state="A"), hsa8794_0)
Initial(hsa8794(bf=None, gene="on", state="I"), hsa8794_0)
Initial(hsa8794(bf=None, gene="off", state="I"), hsa8794_0)
Initial(hsa8795(bf=None, gene="on", state="A"), hsa8795_0)
Initial(hsa8795(bf=None, gene="protein", state="A"), hsa8795_0)
Initial(hsa8795(bf=None, gene="on", state="I"), hsa8795_0)
Initial(hsa8795(bf=None, gene="off", state="I"), hsa8795_0)
Initial(hsa27113(bf=None, gene="on", state="A"), hsa27113_0)
Initial(hsa27113(bf=None, gene="protein", state="A"), hsa27113_0)
Initial(hsa27113(bf=None, gene="on", state="I"), hsa27113_0)
Initial(hsa27113(bf=None, gene="off", state="I"), hsa27113_0)
Initial(hsa2021(bf=None, gene="on", state="A"), hsa2021_0)
Initial(hsa2021(bf=None, gene="protein", state="A"), hsa2021_0)
Initial(hsa2021(bf=None, gene="on", state="I"), hsa2021_0)
Initial(hsa2021(bf=None, gene="off", state="I"), hsa2021_0)
Initial(hsa153090(bf=None, gene="on", state="A"), hsa153090_0)
Initial(hsa153090(bf=None, gene="protein", state="A"), hsa153090_0)
Initial(hsa153090(bf=None, gene="on", state="I"), hsa153090_0)
Initial(hsa153090(bf=None, gene="off", state="I"), hsa153090_0)
Initial(hsa5783(bf=None, gene="on", state="A"), hsa5783_0)
Initial(hsa5783(bf=None, gene="protein", state="A"), hsa5783_0)
Initial(hsa5783(bf=None, gene="on", state="I"), hsa5783_0)
Initial(hsa5783(bf=None, gene="off", state="I"), hsa5783_0)
Initial(hsa7124(bf=None, gene="on", state="A"), hsa7124_0)
Initial(hsa7124(bf=None, gene="protein", state="A"), hsa7124_0)
Initial(hsa7124(bf=None, gene="on", state="I"), hsa7124_0)
Initial(hsa7124(bf=None, gene="off", state="I"), hsa7124_0)
Initial(hsa1509(bf=None, gene="on", state="A"), hsa1509_0)
Initial(hsa1509(bf=None, gene="protein", state="A"), hsa1509_0)
Initial(hsa1509(bf=None, gene="on", state="I"), hsa1509_0)
Initial(hsa1509(bf=None, gene="off", state="I"), hsa1509_0)
Initial(hsa10039(bf=None, gene="on", state="A"), hsa10039_0)
Initial(hsa10039(bf=None, gene="protein", state="A"), hsa10039_0)
Initial(hsa10039(bf=None, gene="on", state="I"), hsa10039_0)
Initial(hsa10039(bf=None, gene="off", state="I"), hsa10039_0)
Initial(hsa10038(bf=None, gene="on", state="A"), hsa10038_0)
Initial(hsa10038(bf=None, gene="protein", state="A"), hsa10038_0)
Initial(hsa10038(bf=None, gene="on", state="I"), hsa10038_0)
Initial(hsa10038(bf=None, gene="off", state="I"), hsa10038_0)
Initial(hsa9131(bf=None, gene="on", state="A"), hsa9131_0)
Initial(hsa9131(bf=None, gene="protein", state="A"), hsa9131_0)
Initial(hsa9131(bf=None, gene="on", state="I"), hsa9131_0)
Initial(hsa9131(bf=None, gene="off", state="I"), hsa9131_0)
Initial(hsa4616(bf=None, gene="on", state="A"), hsa4616_0)
Initial(hsa4616(bf=None, gene="protein", state="A"), hsa4616_0)
Initial(hsa4616(bf=None, gene="on", state="I"), hsa4616_0)
Initial(hsa4616(bf=None, gene="off", state="I"), hsa4616_0)
Initial(hsa8722(bf=None, gene="on", state="A"), hsa8722_0)
Initial(hsa8722(bf=None, gene="protein", state="A"), hsa8722_0)
Initial(hsa8722(bf=None, gene="on", state="I"), hsa8722_0)
Initial(hsa8722(bf=None, gene="off", state="I"), hsa8722_0)
Initial(hsa3845(bf=None, gene="on", state="A"), hsa3845_0)
Initial(hsa3845(bf=None, gene="protein", state="A"), hsa3845_0)
Initial(hsa3845(bf=None, gene="on", state="I"), hsa3845_0)
Initial(hsa3845(bf=None, gene="off", state="I"), hsa3845_0)
Initial(hsa142(bf=None, gene="on", state="A"), hsa142_0)
Initial(hsa142(bf=None, gene="protein", state="A"), hsa142_0)
Initial(hsa142(bf=None, gene="on", state="I"), hsa142_0)
Initial(hsa142(bf=None, gene="off", state="I"), hsa142_0)
Initial(hsa143(bf=None, gene="on", state="A"), hsa143_0)
Initial(hsa143(bf=None, gene="protein", state="A"), hsa143_0)
Initial(hsa143(bf=None, gene="on", state="I"), hsa143_0)
Initial(hsa143(bf=None, gene="off", state="I"), hsa143_0)
Initial(hsa112714(bf=None, gene="on", state="A"), hsa112714_0)
Initial(hsa112714(bf=None, gene="protein", state="A"), hsa112714_0)
Initial(hsa112714(bf=None, gene="on", state="I"), hsa112714_0)
Initial(hsa112714(bf=None, gene="off", state="I"), hsa112714_0)
Initial(hsa10912(bf=None, gene="on", state="A"), hsa10912_0)
Initial(hsa10912(bf=None, gene="protein", state="A"), hsa10912_0)
Initial(hsa10912(bf=None, gene="on", state="I"), hsa10912_0)
Initial(hsa10912(bf=None, gene="off", state="I"), hsa10912_0)
Initial(hsa1513(bf=None, gene="on", state="A"), hsa1513_0)
Initial(hsa1513(bf=None, gene="protein", state="A"), hsa1513_0)
Initial(hsa1513(bf=None, gene="on", state="I"), hsa1513_0)
Initial(hsa1513(bf=None, gene="off", state="I"), hsa1513_0)
Initial(hsa1512(bf=None, gene="on", state="A"), hsa1512_0)
Initial(hsa1512(bf=None, gene="protein", state="A"), hsa1512_0)
Initial(hsa1512(bf=None, gene="on", state="I"), hsa1512_0)
Initial(hsa1512(bf=None, gene="off", state="I"), hsa1512_0)
Initial(hsa1515(bf=None, gene="on", state="A"), hsa1515_0)
Initial(hsa1515(bf=None, gene="protein", state="A"), hsa1515_0)
Initial(hsa1515(bf=None, gene="on", state="I"), hsa1515_0)
Initial(hsa1515(bf=None, gene="off", state="I"), hsa1515_0)
Initial(hsa1514(bf=None, gene="on", state="A"), hsa1514_0)
Initial(hsa1514(bf=None, gene="protein", state="A"), hsa1514_0)
Initial(hsa1514(bf=None, gene="on", state="I"), hsa1514_0)
Initial(hsa1514(bf=None, gene="off", state="I"), hsa1514_0)
Initial(hsa1519(bf=None, gene="on", state="A"), hsa1519_0)
Initial(hsa1519(bf=None, gene="protein", state="A"), hsa1519_0)
Initial(hsa1519(bf=None, gene="on", state="I"), hsa1519_0)
Initial(hsa1519(bf=None, gene="off", state="I"), hsa1519_0)
Initial(hsa84790(bf=None, gene="on", state="A"), hsa84790_0)
Initial(hsa84790(bf=None, gene="protein", state="A"), hsa84790_0)
Initial(hsa84790(bf=None, gene="on", state="I"), hsa84790_0)
Initial(hsa84790(bf=None, gene="off", state="I"), hsa84790_0)
Initial(hsa7132(bf=None, gene="on", state="A"), hsa7132_0)
Initial(hsa7132(bf=None, gene="protein", state="A"), hsa7132_0)
Initial(hsa7132(bf=None, gene="on", state="I"), hsa7132_0)
Initial(hsa7132(bf=None, gene="off", state="I"), hsa7132_0)
Initial(hsa1147(bf=None, gene="on", state="A"), hsa1147_0)
Initial(hsa1147(bf=None, gene="protein", state="A"), hsa1147_0)
Initial(hsa1147(bf=None, gene="on", state="I"), hsa1147_0)
Initial(hsa1147(bf=None, gene="off", state="I"), hsa1147_0)
Initial(hsa8717(bf=None, gene="on", state="A"), hsa8717_0)
Initial(hsa8717(bf=None, gene="protein", state="A"), hsa8717_0)
Initial(hsa8717(bf=None, gene="on", state="I"), hsa8717_0)
Initial(hsa8717(bf=None, gene="off", state="I"), hsa8717_0)
Initial(hsa355(bf=None, gene="on", state="A"), hsa355_0)
Initial(hsa355(bf=None, gene="protein", state="A"), hsa355_0)
Initial(hsa355(bf=None, gene="on", state="I"), hsa355_0)
Initial(hsa355(bf=None, gene="off", state="I"), hsa355_0)
Initial(hsa5894(bf=None, gene="on", state="A"), hsa5894_0)
Initial(hsa5894(bf=None, gene="protein", state="A"), hsa5894_0)
Initial(hsa5894(bf=None, gene="on", state="I"), hsa5894_0)
Initial(hsa5894(bf=None, gene="off", state="I"), hsa5894_0)
Initial(cpdC00076(bf=None), cpdC00076_0)
Initial(hsa9020(bf=None, gene="on", state="A"), hsa9020_0)
Initial(hsa9020(bf=None, gene="protein", state="A"), hsa9020_0)
Initial(hsa9020(bf=None, gene="on", state="I"), hsa9020_0)
Initial(hsa9020(bf=None, gene="off", state="I"), hsa9020_0)
Initial(hsa3551(bf=None, gene="on", state="A"), hsa3551_0)
Initial(hsa3551(bf=None, gene="protein", state="A"), hsa3551_0)
Initial(hsa3551(bf=None, gene="on", state="I"), hsa3551_0)
Initial(hsa3551(bf=None, gene="off", state="I"), hsa3551_0)
Initial(hsa8772(bf=None, gene="on", state="A"), hsa8772_0)
Initial(hsa8772(bf=None, gene="protein", state="A"), hsa8772_0)
Initial(hsa8772(bf=None, gene="on", state="I"), hsa8772_0)
Initial(hsa8772(bf=None, gene="off", state="I"), hsa8772_0)
Initial(hsa572(bf=None, gene="on", state="A"), hsa572_0)
Initial(hsa572(bf=None, gene="protein", state="A"), hsa572_0)
Initial(hsa572(bf=None, gene="on", state="I"), hsa572_0)
Initial(hsa572(bf=None, gene="off", state="I"), hsa572_0)
Initial(hsa356(bf=None, gene="on", state="A"), hsa356_0)
Initial(hsa356(bf=None, gene="protein", state="A"), hsa356_0)
Initial(hsa356(bf=None, gene="on", state="I"), hsa356_0)
Initial(hsa356(bf=None, gene="off", state="I"), hsa356_0)
Initial(hsa4803(bf=None, gene="on", state="A"), hsa4803_0)
Initial(hsa4803(bf=None, gene="protein", state="A"), hsa4803_0)
Initial(hsa4803(bf=None, gene="on", state="I"), hsa4803_0)
Initial(hsa4803(bf=None, gene="off", state="I"), hsa4803_0)
Initial(hsa5366(bf=None, gene="on", state="A"), hsa5366_0)
Initial(hsa5366(bf=None, gene="protein", state="A"), hsa5366_0)
Initial(hsa5366(bf=None, gene="on", state="I"), hsa5366_0)
Initial(hsa5366(bf=None, gene="off", state="I"), hsa5366_0)
Initial(hsa468(bf=None, gene="on", state="A"), hsa468_0)
Initial(hsa468(bf=None, gene="protein", state="A"), hsa468_0)
Initial(hsa468(bf=None, gene="on", state="I"), hsa468_0)
Initial(hsa468(bf=None, gene="off", state="I"), hsa468_0)
Initial(hsa113457(bf=None, gene="on", state="A"), hsa113457_0)
Initial(hsa113457(bf=None, gene="protein", state="A"), hsa113457_0)
Initial(hsa113457(bf=None, gene="on", state="I"), hsa113457_0)
Initial(hsa113457(bf=None, gene="off", state="I"), hsa113457_0)
Initial(hsa27429(bf=None, gene="on", state="A"), hsa27429_0)
Initial(hsa27429(bf=None, gene="protein", state="A"), hsa27429_0)
Initial(hsa27429(bf=None, gene="on", state="I"), hsa27429_0)
Initial(hsa27429(bf=None, gene="off", state="I"), hsa27429_0)
Initial(hsa1522(bf=None, gene="on", state="A"), hsa1522_0)
Initial(hsa1522(bf=None, gene="protein", state="A"), hsa1522_0)
Initial(hsa1522(bf=None, gene="on", state="I"), hsa1522_0)
Initial(hsa1522(bf=None, gene="off", state="I"), hsa1522_0)
Initial(hsa4000(bf=None, gene="on", state="A"), hsa4000_0)
Initial(hsa4000(bf=None, gene="protein", state="A"), hsa4000_0)
Initial(hsa4000(bf=None, gene="on", state="I"), hsa4000_0)
Initial(hsa4000(bf=None, gene="off", state="I"), hsa4000_0)
Initial(hsa4001(bf=None, gene="on", state="A"), hsa4001_0)
Initial(hsa4001(bf=None, gene="protein", state="A"), hsa4001_0)
Initial(hsa4001(bf=None, gene="on", state="I"), hsa4001_0)
Initial(hsa4001(bf=None, gene="off", state="I"), hsa4001_0)
Initial(cpdC00319(bf=None), cpdC00319_0)
Initial(hsa5551(bf=None, gene="on", state="A"), hsa5551_0)
Initial(hsa5551(bf=None, gene="protein", state="A"), hsa5551_0)
Initial(hsa5551(bf=None, gene="on", state="I"), hsa5551_0)
Initial(hsa5551(bf=None, gene="off", state="I"), hsa5551_0)
Initial(hsa2353(bf=None, gene="on", state="A"), hsa2353_0)
Initial(hsa2353(bf=None, gene="protein", state="A"), hsa2353_0)
Initial(hsa2353(bf=None, gene="on", state="I"), hsa2353_0)
Initial(hsa2353(bf=None, gene="off", state="I"), hsa2353_0)
Initial(hsa7186(bf=None, gene="on", state="A"), hsa7186_0)
Initial(hsa7186(bf=None, gene="protein", state="A"), hsa7186_0)
Initial(hsa7186(bf=None, gene="on", state="I"), hsa7186_0)
Initial(hsa7186(bf=None, gene="off", state="I"), hsa7186_0)
Initial(hsa7185(bf=None, gene="on", state="A"), hsa7185_0)
Initial(hsa7185(bf=None, gene="protein", state="A"), hsa7185_0)
Initial(hsa7185(bf=None, gene="on", state="I"), hsa7185_0)
Initial(hsa7185(bf=None, gene="off", state="I"), hsa7185_0)
Initial(hsa8517(bf=None, gene="on", state="A"), hsa8517_0)
Initial(hsa8517(bf=None, gene="protein", state="A"), hsa8517_0)
Initial(hsa8517(bf=None, gene="on", state="I"), hsa8517_0)
Initial(hsa8517(bf=None, gene="off", state="I"), hsa8517_0)
Initial(hsa54205(bf=None, gene="on", state="A"), hsa54205_0)
Initial(hsa54205(bf=None, gene="protein", state="A"), hsa54205_0)
Initial(hsa54205(bf=None, gene="on", state="I"), hsa54205_0)
Initial(hsa54205(bf=None, gene="off", state="I"), hsa54205_0)
Initial(hsa4792(bf=None, gene="on", state="A"), hsa4792_0)
Initial(hsa4792(bf=None, gene="protein", state="A"), hsa4792_0)
Initial(hsa4792(bf=None, gene="on", state="I"), hsa4792_0)
Initial(hsa4792(bf=None, gene="off", state="I"), hsa4792_0)
Initial(hsa4790(bf=None, gene="on", state="A"), hsa4790_0)
Initial(hsa4790(bf=None, gene="protein", state="A"), hsa4790_0)
Initial(hsa4790(bf=None, gene="on", state="I"), hsa4790_0)
Initial(hsa4790(bf=None, gene="off", state="I"), hsa4790_0)
Initial(hsa208(bf=None, gene="on", state="A"), hsa208_0)
Initial(hsa208(bf=None, gene="protein", state="A"), hsa208_0)
Initial(hsa208(bf=None, gene="on", state="I"), hsa208_0)
Initial(hsa208(bf=None, gene="off", state="I"), hsa208_0)
Initial(hsa5970(bf=None, gene="on", state="A"), hsa5970_0)
Initial(hsa5970(bf=None, gene="protein", state="A"), hsa5970_0)
Initial(hsa5970(bf=None, gene="on", state="I"), hsa5970_0)
Initial(hsa5970(bf=None, gene="off", state="I"), hsa5970_0)
Initial(hsa207(bf=None, gene="on", state="A"), hsa207_0)
Initial(hsa207(bf=None, gene="protein", state="A"), hsa207_0)
Initial(hsa207(bf=None, gene="on", state="I"), hsa207_0)
Initial(hsa207(bf=None, gene="off", state="I"), hsa207_0)
Initial(hsa100506742(bf=None, gene="on", state="A"), hsa100506742_0)
Initial(hsa100506742(bf=None, gene="protein", state="A"), hsa100506742_0)
Initial(hsa100506742(bf=None, gene="on", state="I"), hsa100506742_0)
Initial(hsa100506742(bf=None, gene="off", state="I"), hsa100506742_0)
Initial(hsa1647(bf=None, gene="on", state="A"), hsa1647_0)
Initial(hsa1647(bf=None, gene="protein", state="A"), hsa1647_0)
Initial(hsa1647(bf=None, gene="on", state="I"), hsa1647_0)
Initial(hsa1647(bf=None, gene="off", state="I"), hsa1647_0)
Initial(hsa23533(bf=None, gene="on", state="A"), hsa23533_0)
Initial(hsa23533(bf=None, gene="protein", state="A"), hsa23533_0)
Initial(hsa23533(bf=None, gene="on", state="I"), hsa23533_0)
Initial(hsa23533(bf=None, gene="off", state="I"), hsa23533_0)
Initial(hsa1439(bf=None, gene="on", state="A"), hsa1439_0)
Initial(hsa1439(bf=None, gene="protein", state="A"), hsa1439_0)
Initial(hsa1439(bf=None, gene="on", state="I"), hsa1439_0)
Initial(hsa1439(bf=None, gene="off", state="I"), hsa1439_0)
Initial(hsa1649(bf=None, gene="on", state="A"), hsa1649_0)
Initial(hsa1649(bf=None, gene="protein", state="A"), hsa1649_0)
Initial(hsa1649(bf=None, gene="on", state="I"), hsa1649_0)
Initial(hsa1649(bf=None, gene="off", state="I"), hsa1649_0)
Initial(hsa1508(bf=None, gene="on", state="A"), hsa1508_0)
Initial(hsa1508(bf=None, gene="protein", state="A"), hsa1508_0)
Initial(hsa1508(bf=None, gene="on", state="I"), hsa1508_0)
Initial(hsa1508(bf=None, gene="off", state="I"), hsa1508_0)
Initial(hsa8503(bf=None, gene="on", state="A"), hsa8503_0)
Initial(hsa8503(bf=None, gene="protein", state="A"), hsa8503_0)
Initial(hsa8503(bf=None, gene="on", state="I"), hsa8503_0)
Initial(hsa8503(bf=None, gene="off", state="I"), hsa8503_0)
Initial(hsa598(bf=None, gene="on", state="A"), hsa598_0)
Initial(hsa598(bf=None, gene="protein", state="A"), hsa598_0)
Initial(hsa598(bf=None, gene="on", state="I"), hsa598_0)
Initial(hsa598(bf=None, gene="off", state="I"), hsa598_0)
Initial(hsa597(bf=None, gene="on", state="A"), hsa597_0)
Initial(hsa597(bf=None, gene="protein", state="A"), hsa597_0)
Initial(hsa597(bf=None, gene="on", state="I"), hsa597_0)
Initial(hsa597(bf=None, gene="off", state="I"), hsa597_0)
Initial(hsa596(bf=None, gene="on", state="A"), hsa596_0)
Initial(hsa596(bf=None, gene="protein", state="A"), hsa596_0)
Initial(hsa596(bf=None, gene="on", state="I"), hsa596_0)
Initial(hsa596(bf=None, gene="off", state="I"), hsa596_0)
Initial(hsa3562(bf=None, gene="on", state="A"), hsa3562_0)
Initial(hsa3562(bf=None, gene="protein", state="A"), hsa3562_0)
Initial(hsa3562(bf=None, gene="on", state="I"), hsa3562_0)
Initial(hsa3562(bf=None, gene="off", state="I"), hsa3562_0)
Initial(hsa3563(bf=None, gene="on", state="A"), hsa3563_0)
Initial(hsa3563(bf=None, gene="protein", state="A"), hsa3563_0)
Initial(hsa3563(bf=None, gene="on", state="I"), hsa3563_0)
Initial(hsa3563(bf=None, gene="off", state="I"), hsa3563_0)
Initial(hsa841(bf=None, gene="on", state="A"), hsa841_0)
Initial(hsa841(bf=None, gene="protein", state="A"), hsa841_0)
Initial(hsa841(bf=None, gene="on", state="I"), hsa841_0)
Initial(hsa841(bf=None, gene="off", state="I"), hsa841_0)
Initial(hsa840(bf=None, gene="on", state="A"), hsa840_0)
Initial(hsa840(bf=None, gene="protein", state="A"), hsa840_0)
Initial(hsa840(bf=None, gene="on", state="I"), hsa840_0)
Initial(hsa840(bf=None, gene="off", state="I"), hsa840_0)
Initial(hsa843(bf=None, gene="on", state="A"), hsa843_0)
Initial(hsa843(bf=None, gene="protein", state="A"), hsa843_0)
Initial(hsa843(bf=None, gene="on", state="I"), hsa843_0)
Initial(hsa843(bf=None, gene="off", state="I"), hsa843_0)
Initial(hsa842(bf=None, gene="on", state="A"), hsa842_0)
Initial(hsa842(bf=None, gene="protein", state="A"), hsa842_0)
Initial(hsa842(bf=None, gene="on", state="I"), hsa842_0)
Initial(hsa842(bf=None, gene="off", state="I"), hsa842_0)
Initial(hsa1676(bf=None, gene="on", state="A"), hsa1676_0)
Initial(hsa1676(bf=None, gene="protein", state="A"), hsa1676_0)
Initial(hsa1676(bf=None, gene="on", state="I"), hsa1676_0)
Initial(hsa1676(bf=None, gene="off", state="I"), hsa1676_0)
Initial(hsa1677(bf=None, gene="on", state="A"), hsa1677_0)
Initial(hsa1677(bf=None, gene="protein", state="A"), hsa1677_0)
Initial(hsa1677(bf=None, gene="on", state="I"), hsa1677_0)
Initial(hsa1677(bf=None, gene="off", state="I"), hsa1677_0)
Initial(hsa1520(bf=None, gene="on", state="A"), hsa1520_0)
Initial(hsa1520(bf=None, gene="protein", state="A"), hsa1520_0)
Initial(hsa1520(bf=None, gene="on", state="I"), hsa1520_0)
Initial(hsa1520(bf=None, gene="off", state="I"), hsa1520_0)
Initial(hsa4893(bf=None, gene="on", state="A"), hsa4893_0)
Initial(hsa4893(bf=None, gene="protein", state="A"), hsa4893_0)
Initial(hsa4893(bf=None, gene="on", state="I"), hsa4893_0)
Initial(hsa4893(bf=None, gene="off", state="I"), hsa4893_0)
Initial(hsa1075(bf=None, gene="on", state="A"), hsa1075_0)
Initial(hsa1075(bf=None, gene="protein", state="A"), hsa1075_0)
Initial(hsa1075(bf=None, gene="on", state="I"), hsa1075_0)
Initial(hsa1075(bf=None, gene="off", state="I"), hsa1075_0)
Initial(cpdC05981(bf=None), cpdC05981_0)
Initial(hsa472(bf=None, gene="on", state="A"), hsa472_0)
Initial(hsa472(bf=None, gene="protein", state="A"), hsa472_0)
Initial(hsa472(bf=None, gene="on", state="I"), hsa472_0)
Initial(hsa472(bf=None, gene="off", state="I"), hsa472_0)
Initial(hsa63970(bf=None, gene="on", state="A"), hsa63970_0)
Initial(hsa63970(bf=None, gene="protein", state="A"), hsa63970_0)
Initial(hsa63970(bf=None, gene="on", state="I"), hsa63970_0)
Initial(hsa63970(bf=None, gene="off", state="I"), hsa63970_0)
Initial(hsa5170(bf=None, gene="on", state="A"), hsa5170_0)
Initial(hsa5170(bf=None, gene="protein", state="A"), hsa5170_0)
Initial(hsa5170(bf=None, gene="on", state="I"), hsa5170_0)
Initial(hsa5170(bf=None, gene="off", state="I"), hsa5170_0)
Initial(hsa3002(bf=None, gene="on", state="A"), hsa3002_0)
Initial(hsa3002(bf=None, gene="protein", state="A"), hsa3002_0)
Initial(hsa3002(bf=None, gene="on", state="I"), hsa3002_0)
Initial(hsa3002(bf=None, gene="off", state="I"), hsa3002_0)
Initial(hsa5595(bf=None, gene="on", state="A"), hsa5595_0)
Initial(hsa5595(bf=None, gene="protein", state="A"), hsa5595_0)
Initial(hsa5595(bf=None, gene="on", state="I"), hsa5595_0)
Initial(hsa5595(bf=None, gene="off", state="I"), hsa5595_0)
Initial(hsa5594(bf=None, gene="on", state="A"), hsa5594_0)
Initial(hsa5594(bf=None, gene="protein", state="A"), hsa5594_0)
Initial(hsa5594(bf=None, gene="on", state="I"), hsa5594_0)
Initial(hsa5594(bf=None, gene="off", state="I"), hsa5594_0)
Initial(hsa2081(bf=None, gene="on", state="A"), hsa2081_0)
Initial(hsa2081(bf=None, gene="protein", state="A"), hsa2081_0)
Initial(hsa2081(bf=None, gene="on", state="I"), hsa2081_0)
Initial(hsa2081(bf=None, gene="off", state="I"), hsa2081_0)
Initial(hsa5290(bf=None, gene="on", state="A"), hsa5290_0)
Initial(hsa5290(bf=None, gene="protein", state="A"), hsa5290_0)
Initial(hsa5290(bf=None, gene="on", state="I"), hsa5290_0)
Initial(hsa5290(bf=None, gene="off", state="I"), hsa5290_0)
Initial(hsa5291(bf=None, gene="on", state="A"), hsa5291_0)
Initial(hsa5291(bf=None, gene="protein", state="A"), hsa5291_0)
Initial(hsa5291(bf=None, gene="on", state="I"), hsa5291_0)
Initial(hsa5291(bf=None, gene="off", state="I"), hsa5291_0)
Initial(hsa5293(bf=None, gene="on", state="A"), hsa5293_0)
Initial(hsa5293(bf=None, gene="protein", state="A"), hsa5293_0)
Initial(hsa5293(bf=None, gene="on", state="I"), hsa5293_0)
Initial(hsa5293(bf=None, gene="off", state="I"), hsa5293_0)
Initial(hsa5294(bf=None, gene="on", state="A"), hsa5294_0)
Initial(hsa5294(bf=None, gene="protein", state="A"), hsa5294_0)
Initial(hsa5294(bf=None, gene="on", state="I"), hsa5294_0)
Initial(hsa5294(bf=None, gene="off", state="I"), hsa5294_0)
Initial(hsa581(bf=None, gene="on", state="A"), hsa581_0)
Initial(hsa581(bf=None, gene="protein", state="A"), hsa581_0)
Initial(hsa581(bf=None, gene="on", state="I"), hsa581_0)
Initial(hsa581(bf=None, gene="off", state="I"), hsa581_0)
Initial(hsa5296(bf=None, gene="on", state="A"), hsa5296_0)
Initial(hsa5296(bf=None, gene="protein", state="A"), hsa5296_0)
Initial(hsa5296(bf=None, gene="on", state="I"), hsa5296_0)
Initial(hsa5296(bf=None, gene="off", state="I"), hsa5296_0)
Initial(hsa578(bf=None, gene="on", state="A"), hsa578_0)
Initial(hsa578(bf=None, gene="protein", state="A"), hsa578_0)
Initial(hsa578(bf=None, gene="on", state="I"), hsa578_0)
Initial(hsa578(bf=None, gene="off", state="I"), hsa578_0)
Initial(hsa7278(bf=None, gene="on", state="A"), hsa7278_0)
Initial(hsa7278(bf=None, gene="protein", state="A"), hsa7278_0)
Initial(hsa7278(bf=None, gene="on", state="I"), hsa7278_0)
Initial(hsa7278(bf=None, gene="off", state="I"), hsa7278_0)
Initial(hsa4170(bf=None, gene="on", state="A"), hsa4170_0)
Initial(hsa4170(bf=None, gene="protein", state="A"), hsa4170_0)
Initial(hsa4170(bf=None, gene="on", state="I"), hsa4170_0)
Initial(hsa4170(bf=None, gene="off", state="I"), hsa4170_0)
Initial(hsa7277(bf=None, gene="on", state="A"), hsa7277_0)
Initial(hsa7277(bf=None, gene="protein", state="A"), hsa7277_0)
Initial(hsa7277(bf=None, gene="on", state="I"), hsa7277_0)
Initial(hsa7277(bf=None, gene="off", state="I"), hsa7277_0)
Initial(cpdC00027(bf=None), cpdC00027_0)
Initial(hsa7846(bf=None, gene="on", state="A"), hsa7846_0)
Initial(hsa7846(bf=None, gene="protein", state="A"), hsa7846_0)
Initial(hsa7846(bf=None, gene="on", state="I"), hsa7846_0)
Initial(hsa7846(bf=None, gene="off", state="I"), hsa7846_0)
Initial(hsa10018(bf=None, gene="on", state="A"), hsa10018_0)
Initial(hsa10018(bf=None, gene="protein", state="A"), hsa10018_0)
Initial(hsa10018(bf=None, gene="on", state="I"), hsa10018_0)
Initial(hsa10018(bf=None, gene="off", state="I"), hsa10018_0)
Initial(hsa10376(bf=None, gene="on", state="A"), hsa10376_0)
Initial(hsa10376(bf=None, gene="protein", state="A"), hsa10376_0)
Initial(hsa10376(bf=None, gene="on", state="I"), hsa10376_0)
Initial(hsa10376(bf=None, gene="off", state="I"), hsa10376_0)
Initial(hsa55367(bf=None, gene="on", state="A"), hsa55367_0)
Initial(hsa55367(bf=None, gene="protein", state="A"), hsa55367_0)
Initial(hsa55367(bf=None, gene="on", state="I"), hsa55367_0)
Initial(hsa55367(bf=None, gene="off", state="I"), hsa55367_0)
Initial(hsa3708(bf=None, gene="on", state="A"), hsa3708_0)
Initial(hsa3708(bf=None, gene="protein", state="A"), hsa3708_0)
Initial(hsa3708(bf=None, gene="on", state="I"), hsa3708_0)
Initial(hsa3708(bf=None, gene="off", state="I"), hsa3708_0)
Initial(hsa3709(bf=None, gene="on", state="A"), hsa3709_0)
Initial(hsa3709(bf=None, gene="protein", state="A"), hsa3709_0)
Initial(hsa3709(bf=None, gene="on", state="I"), hsa3709_0)
Initial(hsa3709(bf=None, gene="off", state="I"), hsa3709_0)
Initial(hsa10000(bf=None, gene="on", state="A"), hsa10000_0)
Initial(hsa10000(bf=None, gene="protein", state="A"), hsa10000_0)
Initial(hsa10000(bf=None, gene="on", state="I"), hsa10000_0)
Initial(hsa10000(bf=None, gene="off", state="I"), hsa10000_0)
Initial(hsa56616(bf=None, gene="on", state="A"), hsa56616_0)
Initial(hsa56616(bf=None, gene="protein", state="A"), hsa56616_0)
Initial(hsa56616(bf=None, gene="on", state="I"), hsa56616_0)
Initial(hsa56616(bf=None, gene="off", state="I"), hsa56616_0)
Initial(hsa84823(bf=None, gene="on", state="A"), hsa84823_0)
Initial(hsa84823(bf=None, gene="protein", state="A"), hsa84823_0)
Initial(hsa84823(bf=None, gene="on", state="I"), hsa84823_0)
Initial(hsa84823(bf=None, gene="off", state="I"), hsa84823_0)
Initial(hsa6709(bf=None, gene="on", state="A"), hsa6709_0)
Initial(hsa6709(bf=None, gene="protein", state="A"), hsa6709_0)
Initial(hsa6709(bf=None, gene="on", state="I"), hsa6709_0)
Initial(hsa6709(bf=None, gene="off", state="I"), hsa6709_0)
Initial(hsa6708(bf=None, gene="on", state="A"), hsa6708_0)
Initial(hsa6708(bf=None, gene="protein", state="A"), hsa6708_0)
Initial(hsa6708(bf=None, gene="on", state="I"), hsa6708_0)
Initial(hsa6708(bf=None, gene="off", state="I"), hsa6708_0)
Initial(hsa332(bf=None, gene="on", state="A"), hsa332_0)
Initial(hsa332(bf=None, gene="protein", state="A"), hsa332_0)
Initial(hsa332(bf=None, gene="on", state="I"), hsa332_0)
Initial(hsa332(bf=None, gene="off", state="I"), hsa332_0)
Initial(hsa331(bf=None, gene="on", state="A"), hsa331_0)
Initial(hsa331(bf=None, gene="protein", state="A"), hsa331_0)
Initial(hsa331(bf=None, gene="on", state="I"), hsa331_0)
Initial(hsa331(bf=None, gene="off", state="I"), hsa331_0)
Initial(hsa330(bf=None, gene="on", state="A"), hsa330_0)
Initial(hsa330(bf=None, gene="protein", state="A"), hsa330_0)
Initial(hsa330(bf=None, gene="on", state="I"), hsa330_0)
Initial(hsa330(bf=None, gene="off", state="I"), hsa330_0)
Initial(hsa824(bf=None, gene="on", state="A"), hsa824_0)
Initial(hsa824(bf=None, gene="protein", state="A"), hsa824_0)
Initial(hsa824(bf=None, gene="on", state="I"), hsa824_0)
Initial(hsa824(bf=None, gene="off", state="I"), hsa824_0)
Initial(hsa823(bf=None, gene="on", state="A"), hsa823_0)
Initial(hsa823(bf=None, gene="protein", state="A"), hsa823_0)
Initial(hsa823(bf=None, gene="on", state="I"), hsa823_0)
Initial(hsa823(bf=None, gene="off", state="I"), hsa823_0)
Initial(hsa8837(bf=None, gene="on", state="A"), hsa8837_0)
Initial(hsa8837(bf=None, gene="protein", state="A"), hsa8837_0)
Initial(hsa8837(bf=None, gene="on", state="I"), hsa8837_0)
Initial(hsa8837(bf=None, gene="off", state="I"), hsa8837_0)
Initial(hsa1616(bf=None, gene="on", state="A"), hsa1616_0)
Initial(hsa1616(bf=None, gene="protein", state="A"), hsa1616_0)
Initial(hsa1616(bf=None, gene="on", state="I"), hsa1616_0)
Initial(hsa1616(bf=None, gene="off", state="I"), hsa1616_0)
Initial(hsa5605(bf=None, gene="on", state="A"), hsa5605_0)
Initial(hsa5605(bf=None, gene="protein", state="A"), hsa5605_0)
Initial(hsa5605(bf=None, gene="on", state="I"), hsa5605_0)
Initial(hsa5605(bf=None, gene="off", state="I"), hsa5605_0)
Initial(hsa5604(bf=None, gene="on", state="A"), hsa5604_0)
Initial(hsa5604(bf=None, gene="protein", state="A"), hsa5604_0)
Initial(hsa5604(bf=None, gene="on", state="I"), hsa5604_0)
Initial(hsa5604(bf=None, gene="off", state="I"), hsa5604_0)
Initial(hsa5602(bf=None, gene="on", state="A"), hsa5602_0)
Initial(hsa5602(bf=None, gene="protein", state="A"), hsa5602_0)
Initial(hsa5602(bf=None, gene="on", state="I"), hsa5602_0)
Initial(hsa5602(bf=None, gene="off", state="I"), hsa5602_0)
Initial(hsa5601(bf=None, gene="on", state="A"), hsa5601_0)
Initial(hsa5601(bf=None, gene="protein", state="A"), hsa5601_0)
Initial(hsa5601(bf=None, gene="on", state="I"), hsa5601_0)
Initial(hsa5601(bf=None, gene="off", state="I"), hsa5601_0)
Initial(hsa8743(bf=None, gene="on", state="A"), hsa8743_0)
Initial(hsa8743(bf=None, gene="protein", state="A"), hsa8743_0)
Initial(hsa8743(bf=None, gene="on", state="I"), hsa8743_0)
Initial(hsa8743(bf=None, gene="off", state="I"), hsa8743_0)
Initial(hsa60(bf=None, gene="on", state="A"), hsa60_0)
Initial(hsa60(bf=None, gene="protein", state="A"), hsa60_0)
Initial(hsa60(bf=None, gene="on", state="I"), hsa60_0)
Initial(hsa60(bf=None, gene="off", state="I"), hsa60_0)
Initial(hsa4914(bf=None, gene="on", state="A"), hsa4914_0)
Initial(hsa4914(bf=None, gene="protein", state="A"), hsa4914_0)
Initial(hsa4914(bf=None, gene="on", state="I"), hsa4914_0)
Initial(hsa4914(bf=None, gene="off", state="I"), hsa4914_0)
Initial(hsa637(bf=None, gene="on", state="A"), hsa637_0)
Initial(hsa637(bf=None, gene="protein", state="A"), hsa637_0)
Initial(hsa637(bf=None, gene="on", state="I"), hsa637_0)
Initial(hsa637(bf=None, gene="off", state="I"), hsa637_0)
Initial(hsa9451(bf=None, gene="on", state="A"), hsa9451_0)
Initial(hsa9451(bf=None, gene="protein", state="A"), hsa9451_0)
Initial(hsa9451(bf=None, gene="on", state="I"), hsa9451_0)
Initial(hsa9451(bf=None, gene="off", state="I"), hsa9451_0)
Initial(hsa3710(bf=None, gene="on", state="A"), hsa3710_0)
Initial(hsa3710(bf=None, gene="protein", state="A"), hsa3710_0)
Initial(hsa3710(bf=None, gene="on", state="I"), hsa3710_0)
Initial(hsa3710(bf=None, gene="off", state="I"), hsa3710_0)
Initial(hsa79861(bf=None, gene="on", state="A"), hsa79861_0)
Initial(hsa79861(bf=None, gene="protein", state="A"), hsa79861_0)
Initial(hsa79861(bf=None, gene="on", state="I"), hsa79861_0)
Initial(hsa79861(bf=None, gene="off", state="I"), hsa79861_0)
Initial(hsa5295(bf=None, gene="on", state="A"), hsa5295_0)
Initial(hsa5295(bf=None, gene="protein", state="A"), hsa5295_0)
Initial(hsa5295(bf=None, gene="on", state="I"), hsa5295_0)
Initial(hsa5295(bf=None, gene="off", state="I"), hsa5295_0)
Initial(hsa329(bf=None, gene="on", state="A"), hsa329_0)
Initial(hsa329(bf=None, gene="protein", state="A"), hsa329_0)
Initial(hsa329(bf=None, gene="on", state="I"), hsa329_0)
Initial(hsa329(bf=None, gene="off", state="I"), hsa329_0)
Initial(hsa5599(bf=None, gene="on", state="A"), hsa5599_0)
Initial(hsa5599(bf=None, gene="protein", state="A"), hsa5599_0)
Initial(hsa5599(bf=None, gene="on", state="I"), hsa5599_0)
Initial(hsa5599(bf=None, gene="off", state="I"), hsa5599_0)
Initial(hsa835(bf=None, gene="on", state="A"), hsa835_0)
Initial(hsa835(bf=None, gene="protein", state="A"), hsa835_0)
Initial(hsa835(bf=None, gene="on", state="I"), hsa835_0)
Initial(hsa835(bf=None, gene="off", state="I"), hsa835_0)
Initial(hsa836(bf=None, gene="on", state="A"), hsa836_0)
Initial(hsa836(bf=None, gene="protein", state="A"), hsa836_0)
Initial(hsa836(bf=None, gene="on", state="I"), hsa836_0)
Initial(hsa836(bf=None, gene="off", state="I"), hsa836_0)
Initial(hsa5414(bf=None, gene="on", state="A"), hsa5414_0)
Initial(hsa5414(bf=None, gene="protein", state="A"), hsa5414_0)
Initial(hsa5414(bf=None, gene="on", state="I"), hsa5414_0)
Initial(hsa5414(bf=None, gene="off", state="I"), hsa5414_0)
Initial(hsa1521(bf=None, gene="on", state="A"), hsa1521_0)
Initial(hsa1521(bf=None, gene="protein", state="A"), hsa1521_0)
Initial(hsa1521(bf=None, gene="on", state="I"), hsa1521_0)
Initial(hsa1521(bf=None, gene="off", state="I"), hsa1521_0)
Initial(hsa839(bf=None, gene="on", state="A"), hsa839_0)
Initial(hsa839(bf=None, gene="protein", state="A"), hsa839_0)
Initial(hsa839(bf=None, gene="on", state="I"), hsa839_0)
Initial(hsa839(bf=None, gene="off", state="I"), hsa839_0)
Initial(hsa7157(bf=None, gene="on", state="A"), hsa7157_0)
Initial(hsa7157(bf=None, gene="protein", state="A"), hsa7157_0)
Initial(hsa7157(bf=None, gene="on", state="I"), hsa7157_0)
Initial(hsa7157(bf=None, gene="off", state="I"), hsa7157_0)
Initial(hsa4217(bf=None, gene="on", state="A"), hsa4217_0)
Initial(hsa4217(bf=None, gene="protein", state="A"), hsa4217_0)
Initial(hsa4217(bf=None, gene="on", state="I"), hsa4217_0)
Initial(hsa4217(bf=None, gene="off", state="I"), hsa4217_0)
Initial(hsa1965(bf=None, gene="on", state="A"), hsa1965_0)
Initial(hsa1965(bf=None, gene="protein", state="A"), hsa1965_0)
Initial(hsa1965(bf=None, gene="on", state="I"), hsa1965_0)
Initial(hsa1965(bf=None, gene="off", state="I"), hsa1965_0)
Initial(hsa8739(bf=None, gene="on", state="A"), hsa8739_0)
Initial(hsa8739(bf=None, gene="protein", state="A"), hsa8739_0)
Initial(hsa8739(bf=None, gene="on", state="I"), hsa8739_0)
Initial(hsa8739(bf=None, gene="off", state="I"), hsa8739_0)
Initial(hsa3265(bf=None, gene="on", state="A"), hsa3265_0)
Initial(hsa3265(bf=None, gene="protein", state="A"), hsa3265_0)
Initial(hsa3265(bf=None, gene="on", state="I"), hsa3265_0)
Initial(hsa3265(bf=None, gene="off", state="I"), hsa3265_0)
Initial(hsa8737(bf=None, gene="on", state="A"), hsa8737_0)
Initial(hsa8737(bf=None, gene="protein", state="A"), hsa8737_0)
Initial(hsa8737(bf=None, gene="on", state="I"), hsa8737_0)
Initial(hsa8737(bf=None, gene="off", state="I"), hsa8737_0)
Initial(hsa71(bf=None, gene="on", state="A"), hsa71_0)
Initial(hsa71(bf=None, gene="protein", state="A"), hsa71_0)
Initial(hsa71(bf=None, gene="on", state="I"), hsa71_0)
Initial(hsa71(bf=None, gene="off", state="I"), hsa71_0)
catalyze(hsa317(gene="protein", state="A"), hsa842(gene="protein", state="I"),
         hsa842(gene="protein", state="A"), [kf0, kr0, kc0])
catalyze(hsa8793(gene="protein", state="A"),
         hsa8772(gene="protein", state="I"),
         hsa8772(gene="protein", state="A"), [kf1, kr1, kc1])
catalyze(hsa3725(gene="protein"), hsa8739(gene="off"),
         hsa8739(gene="on", state="A"), [kf2, kr2, kc2])
Rule("hsa8739_expression_2",
     hsa8739(bf=None, gene="on", state="A") >> hsa8739(bf=None, gene="on",
                                                       state="A") + hsa8739(
         bf=None, gene="protein", state="A"), kr2)
catalyze(hsa3725(gene="protein"), hsa355(gene="off"),
         hsa355(gene="on", state="A"), [kf3, kr3, kc3])
Rule("hsa355_expression_3",
     hsa355(bf=None, gene="on", state="A") >> hsa355(bf=None, gene="on",
                                                     state="A") + hsa355(
         bf=None, gene="protein", state="A"), kr3)
catalyze(hsa3725(gene="protein"), hsa10018(gene="off"),
         hsa10018(gene="on", state="A"), [kf4, kr4, kc4])
Rule("hsa10018_expression_4",
     hsa10018(bf=None, gene="on", state="A") >> hsa10018(bf=None, gene="on",
                                                         state="A") + hsa10018(
         bf=None, gene="protein", state="A"), kr4)
catalyze(hsa3725(gene="protein"), hsa7157(gene="off"),
         hsa7157(gene="on", state="A"), [kf5, kr5, kc5])
Rule("hsa7157_expression_5",
     hsa7157(bf=None, gene="on", state="A") >> hsa7157(bf=None, gene="on",
                                                       state="A") + hsa7157(
         bf=None, gene="protein", state="A"), kr5)
catalyze(hsa3725(gene="protein"), hsa356(gene="off"),
         hsa356(gene="on", state="A"), [kf6, kr6, kc6])
Rule("hsa356_expression_6",
     hsa356(bf=None, gene="on", state="A") >> hsa356(bf=None, gene="on",
                                                     state="A") + hsa356(
         bf=None, gene="protein", state="A"), kr6)
catalyze(hsa8797(gene="protein", state="A"),
         hsa8772(gene="protein", state="I"),
         hsa8772(gene="protein", state="A"), [kf7, kr7, kc7])
catalyze(hsa8794(gene="protein", state="A"),
         hsa8772(gene="protein", state="I"),
         hsa8772(gene="protein", state="A"), [kf8, kr8, kc8])
catalyze(hsa8795(gene="protein", state="A"),
         hsa8772(gene="protein", state="I"),
         hsa8772(gene="protein", state="A"), [kf9, kr9, kc9])
catalyze(hsa27113(gene="protein", state="A"),
         hsa596(gene="protein", state="A"), hsa596(gene="protein", state="I"),
         [kf10, kr10, kc10])
catalyze(hsa27113(gene="protein", state="A"),
         hsa598(gene="protein", state="A"), hsa598(gene="protein", state="I"),
         [kf11, kr11, kc11])
catalyze(hsa153090(gene="protein", state="A"),
         hsa4217(gene="protein", state="I"),
         hsa4217(gene="protein", state="A"), [kf12, kr12, kc12])
catalyze(hsa7124(gene="protein", state="A"),
         hsa7132(gene="protein", state="I"),
         hsa7132(gene="protein", state="A"), [kf13, kr13, kc13])
catalyze(hsa1509(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf14, kr14, kc14])
catalyze(hsa1509(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf15, kr15, kc15])
catalyze(hsa1509(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf16, kr16, kc16])
catalyze(hsa1509(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf17, kr17, kc17])
catalyze(hsa1509(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf18, kr18, kc18])
catalyze(hsa1509(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf19, kr19, kc19])
catalyze(hsa1509(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf20, kr20, kc20])
catalyze(hsa8722(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf21, kr21, kc21])
catalyze(hsa8722(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf22, kr22, kc22])
catalyze(hsa8722(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf23, kr23, kc23])
catalyze(hsa8722(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf24, kr24, kc24])
catalyze(hsa8722(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf25, kr25, kc25])
catalyze(hsa8722(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf26, kr26, kc26])
catalyze(hsa8722(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf27, kr27, kc27])
catalyze(hsa3845(gene="protein", state="A"),
         hsa5894(gene="protein", state="I"),
         hsa5894(gene="protein", state="A"), [kf28, kr28, kc28])
catalyze(hsa1513(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf29, kr29, kc29])
catalyze(hsa1513(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf30, kr30, kc30])
catalyze(hsa1513(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf31, kr31, kc31])
catalyze(hsa1513(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf32, kr32, kc32])
catalyze(hsa1513(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf33, kr33, kc33])
catalyze(hsa1513(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf34, kr34, kc34])
catalyze(hsa1513(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf35, kr35, kc35])
catalyze(hsa1512(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf36, kr36, kc36])
catalyze(hsa1512(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf37, kr37, kc37])
catalyze(hsa1512(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf38, kr38, kc38])
catalyze(hsa1512(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf39, kr39, kc39])
catalyze(hsa1512(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf40, kr40, kc40])
catalyze(hsa1512(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf41, kr41, kc41])
catalyze(hsa1512(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf42, kr42, kc42])
catalyze(hsa1515(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf43, kr43, kc43])
catalyze(hsa1515(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf44, kr44, kc44])
catalyze(hsa1515(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf45, kr45, kc45])
catalyze(hsa1515(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf46, kr46, kc46])
catalyze(hsa1515(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf47, kr47, kc47])
catalyze(hsa1515(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf48, kr48, kc48])
catalyze(hsa1515(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf49, kr49, kc49])
catalyze(hsa1514(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf50, kr50, kc50])
catalyze(hsa1514(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf51, kr51, kc51])
catalyze(hsa1514(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf52, kr52, kc52])
catalyze(hsa1514(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf53, kr53, kc53])
catalyze(hsa1514(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf54, kr54, kc54])
catalyze(hsa1514(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf55, kr55, kc55])
catalyze(hsa1514(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf56, kr56, kc56])
catalyze(hsa1519(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf57, kr57, kc57])
catalyze(hsa1519(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf58, kr58, kc58])
catalyze(hsa1519(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf59, kr59, kc59])
catalyze(hsa1519(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf60, kr60, kc60])
catalyze(hsa1519(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf61, kr61, kc61])
catalyze(hsa1519(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf62, kr62, kc62])
catalyze(hsa1519(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf63, kr63, kc63])
catalyze(hsa1147(gene="protein", state="A"),
         hsa4792(gene="protein", state="I"),
         hsa4792(gene="protein", state="A"), [kf64, kr64, kc64])
catalyze(hsa355(gene="protein", state="A"), hsa8772(gene="protein", state="I"),
         hsa8772(gene="protein", state="A"), [kf65, kr65, kc65])
catalyze(hsa355(gene="protein", state="A"), hsa1616(gene="protein", state="I"),
         hsa1616(gene="protein", state="A"), [kf66, kr66, kc66])
catalyze(hsa5894(gene="protein", state="A"),
         hsa5605(gene="protein", state="I"),
         hsa5605(gene="protein", state="A"), [kf67, kr67, kc67])
catalyze(hsa5894(gene="protein", state="A"),
         hsa5604(gene="protein", state="I"),
         hsa5604(gene="protein", state="A"), [kf68, kr68, kc68])
catalyze(cpdC00076(bf=None), hsa824(gene="protein", state="I"),
         hsa824(gene="protein", state="A"), [kf69, kr69, kc69])
catalyze(cpdC00076(bf=None), hsa823(gene="protein", state="I"),
         hsa823(gene="protein", state="A"), [kf70, kr70, kc70])
catalyze(hsa9020(gene="protein", state="A"),
         hsa3551(gene="protein", state="I"),
         hsa3551(gene="protein", state="A"), [kf71, kr71, kc71])
catalyze(hsa9020(gene="protein", state="A"),
         hsa8517(gene="protein", state="I"),
         hsa8517(gene="protein", state="A"), [kf72, kr72, kc72])
catalyze(hsa9020(gene="protein", state="A"),
         hsa1147(gene="protein", state="I"),
         hsa1147(gene="protein", state="A"), [kf73, kr73, kc73])
catalyze(hsa3551(gene="protein", state="A"),
         hsa4792(gene="protein", state="I"),
         hsa4792(gene="protein", state="A"), [kf74, kr74, kc74])
catalyze(hsa8772(gene="protein", state="A"), hsa841(gene="protein", state="I"),
         hsa841(gene="protein", state="A"), [kf75, kr75, kc75])
catalyze(hsa8772(gene="protein", state="A"), hsa843(gene="protein", state="I"),
         hsa843(gene="protein", state="A"), [kf76, kr76, kc76])
catalyze(hsa572(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf77, kr77, kc77])
catalyze(hsa356(gene="protein", state="A"), hsa355(gene="protein", state="I"),
         hsa355(gene="protein", state="A"), [kf78, kr78, kc78])
catalyze(hsa4803(gene="protein", state="A"),
         hsa4914(gene="protein", state="I"),
         hsa4914(gene="protein", state="A"), [kf79, kr79, kc79])
catalyze(hsa5366(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf80, kr80, kc80])
catalyze(hsa5366(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf81, kr81, kc81])
catalyze(hsa468(gene="protein"), hsa1649(gene="off"),
         hsa1649(gene="on", state="A"), [kf82, kr82, kc82])
Rule("hsa1649_expression_82",
     hsa1649(bf=None, gene="on", state="A") >> hsa1649(bf=None, gene="on",
                                                       state="A") + hsa1649(
         bf=None, gene="protein", state="A"), kr82)
catalyze(hsa27429(gene="protein", state="A"),
         hsa329(gene="protein", state="A"), hsa329(gene="protein", state="I"),
         [kf83, kr83, kc83])
catalyze(hsa27429(gene="protein", state="A"),
         hsa331(gene="protein", state="A"), hsa331(gene="protein", state="I"),
         [kf84, kr84, kc84])
catalyze(hsa27429(gene="protein", state="A"),
         hsa330(gene="protein", state="A"), hsa330(gene="protein", state="I"),
         [kf85, kr85, kc85])
catalyze(hsa27429(gene="protein", state="A"),
         hsa332(gene="protein", state="A"), hsa332(gene="protein", state="I"),
         [kf86, kr86, kc86])
catalyze(hsa1522(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf87, kr87, kc87])
catalyze(hsa1522(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf88, kr88, kc88])
catalyze(hsa1522(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf89, kr89, kc89])
catalyze(hsa1522(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf90, kr90, kc90])
catalyze(hsa1522(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf91, kr91, kc91])
catalyze(hsa1522(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf92, kr92, kc92])
catalyze(hsa1522(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf93, kr93, kc93])
catalyze(hsa2353(gene="protein"), hsa8739(gene="off"),
         hsa8739(gene="on", state="A"), [kf94, kr94, kc94])
Rule("hsa8739_expression_94",
     hsa8739(bf=None, gene="on", state="A") >> hsa8739(bf=None, gene="on",
                                                       state="A") + hsa8739(
         bf=None, gene="protein", state="A"), kr94)
catalyze(hsa2353(gene="protein"), hsa355(gene="off"),
         hsa355(gene="on", state="A"), [kf95, kr95, kc95])
Rule("hsa355_expression_95",
     hsa355(bf=None, gene="on", state="A") >> hsa355(bf=None, gene="on",
                                                     state="A") + hsa355(
         bf=None, gene="protein", state="A"), kr95)
catalyze(hsa2353(gene="protein"), hsa10018(gene="off"),
         hsa10018(gene="on", state="A"), [kf96, kr96, kc96])
Rule("hsa10018_expression_96",
     hsa10018(bf=None, gene="on", state="A") >> hsa10018(bf=None, gene="on",
                                                         state="A") + hsa10018(
         bf=None, gene="protein", state="A"), kr96)
catalyze(hsa2353(gene="protein"), hsa7157(gene="off"),
         hsa7157(gene="on", state="A"), [kf97, kr97, kc97])
Rule("hsa7157_expression_97",
     hsa7157(bf=None, gene="on", state="A") >> hsa7157(bf=None, gene="on",
                                                       state="A") + hsa7157(
         bf=None, gene="protein", state="A"), kr97)
catalyze(hsa2353(gene="protein"), hsa356(gene="off"),
         hsa356(gene="on", state="A"), [kf98, kr98, kc98])
Rule("hsa356_expression_98",
     hsa356(bf=None, gene="on", state="A") >> hsa356(bf=None, gene="on",
                                                     state="A") + hsa356(
         bf=None, gene="protein", state="A"), kr98)
catalyze(hsa7186(gene="protein", state="A"),
         hsa153090(gene="protein", state="I"),
         hsa153090(gene="protein", state="A"), [kf99, kr99, kc99])
catalyze(hsa7186(gene="protein", state="A"),
         hsa9020(gene="protein", state="I"),
         hsa9020(gene="protein", state="A"), [kf100, kr100, kc100])
catalyze(hsa8517(gene="protein", state="A"),
         hsa4792(gene="protein", state="I"),
         hsa4792(gene="protein", state="A"), [kf101, kr101, kc101])
Rule("hsa54205_binds_hsa317_102",
     hsa54205(bf=None, gene="protein") + hsa317(bf=None,
                                                gene="protein") <> hsa54205(
         bf=1, gene="protein") % hsa317(bf=1, gene="protein"), kf102, kr102)
catalyze(hsa54205(gene="protein", state="A"),
         hsa842(gene="protein", state="I"), hsa842(gene="protein", state="A"),
         [kf103, kr103, kc103])
catalyze(hsa4790(gene="protein"), hsa1647(gene="off"),
         hsa1647(gene="on", state="A"), [kf106, kr106, kc106])
Rule("hsa1647_expression_106",
     hsa1647(bf=None, gene="on", state="A") >> hsa1647(bf=None, gene="on",
                                                       state="A") + hsa1647(
         bf=None, gene="protein", state="A"), kr106)
catalyze(hsa4790(gene="protein"), hsa5783(gene="off"),
         hsa5783(gene="on", state="A"), [kf107, kr107, kc107])
Rule("hsa5783_expression_107",
     hsa5783(bf=None, gene="on", state="A") >> hsa5783(bf=None, gene="on",
                                                       state="A") + hsa5783(
         bf=None, gene="protein", state="A"), kr107)
catalyze(hsa4790(gene="protein"), hsa598(gene="off"),
         hsa598(gene="on", state="A"), [kf108, kr108, kc108])
Rule("hsa598_expression_108",
     hsa598(bf=None, gene="on", state="A") >> hsa598(bf=None, gene="on",
                                                     state="A") + hsa598(
         bf=None, gene="protein", state="A"), kr108)
catalyze(hsa4790(gene="protein"), hsa8837(gene="off"),
         hsa8837(gene="on", state="A"), [kf109, kr109, kc109])
Rule("hsa8837_expression_109",
     hsa8837(bf=None, gene="on", state="A") >> hsa8837(bf=None, gene="on",
                                                       state="A") + hsa8837(
         bf=None, gene="protein", state="A"), kr109)
catalyze(hsa4790(gene="protein"), hsa597(gene="off"),
         hsa597(gene="on", state="A"), [kf110, kr110, kc110])
Rule("hsa597_expression_110",
     hsa597(bf=None, gene="on", state="A") >> hsa597(bf=None, gene="on",
                                                     state="A") + hsa597(
         bf=None, gene="protein", state="A"), kr110)
catalyze(hsa4790(gene="protein"), hsa329(gene="off"),
         hsa329(gene="on", state="A"), [kf111, kr111, kc111])
Rule("hsa329_expression_111",
     hsa329(bf=None, gene="on", state="A") >> hsa329(bf=None, gene="on",
                                                     state="A") + hsa329(
         bf=None, gene="protein", state="A"), kr111)
catalyze(hsa4790(gene="protein"), hsa4616(gene="off"),
         hsa4616(gene="on", state="A"), [kf112, kr112, kc112])
Rule("hsa4616_expression_112",
     hsa4616(bf=None, gene="on", state="A") >> hsa4616(bf=None, gene="on",
                                                       state="A") + hsa4616(
         bf=None, gene="protein", state="A"), kr112)
catalyze(hsa4790(gene="protein"), hsa332(gene="off"),
         hsa332(gene="on", state="A"), [kf113, kr113, kc113])
Rule("hsa332_expression_113",
     hsa332(bf=None, gene="on", state="A") >> hsa332(bf=None, gene="on",
                                                     state="A") + hsa332(
         bf=None, gene="protein", state="A"), kr113)
catalyze(hsa4790(gene="protein"), hsa331(gene="off"),
         hsa331(gene="on", state="A"), [kf114, kr114, kc114])
Rule("hsa331_expression_114",
     hsa331(bf=None, gene="on", state="A") >> hsa331(bf=None, gene="on",
                                                     state="A") + hsa331(
         bf=None, gene="protein", state="A"), kr114)
catalyze(hsa4790(gene="protein"), hsa330(gene="off"),
         hsa330(gene="on", state="A"), [kf115, kr115, kc115])
Rule("hsa330_expression_115",
     hsa330(bf=None, gene="on", state="A") >> hsa330(bf=None, gene="on",
                                                     state="A") + hsa330(
         bf=None, gene="protein", state="A"), kr115)
catalyze(hsa4790(gene="protein"), hsa7186(gene="off"),
         hsa7186(gene="on", state="A"), [kf116, kr116, kc116])
Rule("hsa7186_expression_116",
     hsa7186(bf=None, gene="on", state="A") >> hsa7186(bf=None, gene="on",
                                                       state="A") + hsa7186(
         bf=None, gene="protein", state="A"), kr116)
catalyze(hsa4790(gene="protein"), hsa10912(gene="off"),
         hsa10912(gene="on", state="A"), [kf117, kr117, kc117])
Rule("hsa10912_expression_117",
     hsa10912(bf=None, gene="on", state="A") >> hsa10912(bf=None, gene="on",
                                                         state="A") + hsa10912(
         bf=None, gene="protein", state="A"), kr117)
catalyze(hsa4790(gene="protein"), hsa7185(gene="off"),
         hsa7185(gene="on", state="A"), [kf118, kr118, kc118])
Rule("hsa7185_expression_118",
     hsa7185(bf=None, gene="on", state="A") >> hsa7185(bf=None, gene="on",
                                                       state="A") + hsa7185(
         bf=None, gene="protein", state="A"), kr118)
catalyze(hsa208(gene="protein", state="A"), hsa3551(gene="protein", state="I"),
         hsa3551(gene="protein", state="A"), [kf119, kr119, kc119])
catalyze(hsa208(gene="protein", state="A"), hsa1147(gene="protein", state="I"),
         hsa1147(gene="protein", state="A"), [kf120, kr120, kc120])
catalyze(hsa208(gene="protein", state="A"), hsa8517(gene="protein", state="I"),
         hsa8517(gene="protein", state="A"), [kf121, kr121, kc121])
catalyze(hsa208(gene="protein", state="A"), hsa572(gene="protein", state="A"),
         hsa572(gene="protein", state="I"), [kf122, kr122, kc122])
catalyze(hsa5970(gene="protein"), hsa1647(gene="off"),
         hsa1647(gene="on", state="A"), [kf123, kr123, kc123])
Rule("hsa1647_expression_123",
     hsa1647(bf=None, gene="on", state="A") >> hsa1647(bf=None, gene="on",
                                                       state="A") + hsa1647(
         bf=None, gene="protein", state="A"), kr123)
catalyze(hsa5970(gene="protein"), hsa5783(gene="off"),
         hsa5783(gene="on", state="A"), [kf124, kr124, kc124])
Rule("hsa5783_expression_124",
     hsa5783(bf=None, gene="on", state="A") >> hsa5783(bf=None, gene="on",
                                                       state="A") + hsa5783(
         bf=None, gene="protein", state="A"), kr124)
catalyze(hsa5970(gene="protein"), hsa598(gene="off"),
         hsa598(gene="on", state="A"), [kf125, kr125, kc125])
Rule("hsa598_expression_125",
     hsa598(bf=None, gene="on", state="A") >> hsa598(bf=None, gene="on",
                                                     state="A") + hsa598(
         bf=None, gene="protein", state="A"), kr125)
catalyze(hsa5970(gene="protein"), hsa8837(gene="off"),
         hsa8837(gene="on", state="A"), [kf126, kr126, kc126])
Rule("hsa8837_expression_126",
     hsa8837(bf=None, gene="on", state="A") >> hsa8837(bf=None, gene="on",
                                                       state="A") + hsa8837(
         bf=None, gene="protein", state="A"), kr126)
catalyze(hsa5970(gene="protein"), hsa597(gene="off"),
         hsa597(gene="on", state="A"), [kf127, kr127, kc127])
Rule("hsa597_expression_127",
     hsa597(bf=None, gene="on", state="A") >> hsa597(bf=None, gene="on",
                                                     state="A") + hsa597(
         bf=None, gene="protein", state="A"), kr127)
catalyze(hsa5970(gene="protein"), hsa329(gene="off"),
         hsa329(gene="on", state="A"), [kf128, kr128, kc128])
Rule("hsa329_expression_128",
     hsa329(bf=None, gene="on", state="A") >> hsa329(bf=None, gene="on",
                                                     state="A") + hsa329(
         bf=None, gene="protein", state="A"), kr128)
catalyze(hsa5970(gene="protein"), hsa4616(gene="off"),
         hsa4616(gene="on", state="A"), [kf129, kr129, kc129])
Rule("hsa4616_expression_129",
     hsa4616(bf=None, gene="on", state="A") >> hsa4616(bf=None, gene="on",
                                                       state="A") + hsa4616(
         bf=None, gene="protein", state="A"), kr129)
catalyze(hsa5970(gene="protein"), hsa332(gene="off"),
         hsa332(gene="on", state="A"), [kf130, kr130, kc130])
Rule("hsa332_expression_130",
     hsa332(bf=None, gene="on", state="A") >> hsa332(bf=None, gene="on",
                                                     state="A") + hsa332(
         bf=None, gene="protein", state="A"), kr130)
catalyze(hsa5970(gene="protein"), hsa331(gene="off"),
         hsa331(gene="on", state="A"), [kf131, kr131, kc131])
Rule("hsa331_expression_131",
     hsa331(bf=None, gene="on", state="A") >> hsa331(bf=None, gene="on",
                                                     state="A") + hsa331(
         bf=None, gene="protein", state="A"), kr131)
catalyze(hsa5970(gene="protein"), hsa330(gene="off"),
         hsa330(gene="on", state="A"), [kf132, kr132, kc132])
Rule("hsa330_expression_132",
     hsa330(bf=None, gene="on", state="A") >> hsa330(bf=None, gene="on",
                                                     state="A") + hsa330(
         bf=None, gene="protein", state="A"), kr132)
catalyze(hsa5970(gene="protein"), hsa7186(gene="off"),
         hsa7186(gene="on", state="A"), [kf133, kr133, kc133])
Rule("hsa7186_expression_133",
     hsa7186(bf=None, gene="on", state="A") >> hsa7186(bf=None, gene="on",
                                                       state="A") + hsa7186(
         bf=None, gene="protein", state="A"), kr133)
catalyze(hsa5970(gene="protein"), hsa10912(gene="off"),
         hsa10912(gene="on", state="A"), [kf134, kr134, kc134])
Rule("hsa10912_expression_134",
     hsa10912(bf=None, gene="on", state="A") >> hsa10912(bf=None, gene="on",
                                                         state="A") + hsa10912(
         bf=None, gene="protein", state="A"), kr134)
catalyze(hsa5970(gene="protein"), hsa7185(gene="off"),
         hsa7185(gene="on", state="A"), [kf135, kr135, kc135])
Rule("hsa7185_expression_135",
     hsa7185(bf=None, gene="on", state="A") >> hsa7185(bf=None, gene="on",
                                                       state="A") + hsa7185(
         bf=None, gene="protein", state="A"), kr135)
catalyze(hsa207(gene="protein", state="A"), hsa3551(gene="protein", state="I"),
         hsa3551(gene="protein", state="A"), [kf136, kr136, kc136])
catalyze(hsa207(gene="protein", state="A"), hsa1147(gene="protein", state="I"),
         hsa1147(gene="protein", state="A"), [kf137, kr137, kc137])
catalyze(hsa207(gene="protein", state="A"), hsa8517(gene="protein", state="I"),
         hsa8517(gene="protein", state="A"), [kf138, kr138, kc138])
catalyze(hsa207(gene="protein", state="A"), hsa572(gene="protein", state="A"),
         hsa572(gene="protein", state="I"), [kf139, kr139, kc139])
catalyze(hsa100506742(gene="protein", state="A"),
         hsa836(gene="protein", state="I"), hsa836(gene="protein", state="A"),
         [kf140, kr140, kc140])
catalyze(hsa100506742(gene="protein", state="A"),
         hsa840(gene="protein", state="I"), hsa840(gene="protein", state="A"),
         [kf141, kr141, kc141])
catalyze(cpdC05981(bf=None), hsa23533(gene="protein", state="I"),
         hsa23533(gene="protein", state="A"), [kf142, kr142, kc142])
catalyze(hsa1439(gene="protein", state="A"),
         hsa23533(gene="protein", state="I"),
         hsa23533(gene="protein", state="A"), [kf143, kr143, kc143])
catalyze(hsa1439(gene="protein", state="A"),
         hsa8503(gene="protein", state="I"),
         hsa8503(gene="protein", state="A"), [kf144, kr144, kc144])
catalyze(hsa1439(gene="protein", state="A"),
         hsa3265(gene="protein", state="I"),
         hsa3265(gene="protein", state="A"), [kf145, kr145, kc145])
catalyze(hsa1439(gene="protein", state="A"),
         hsa5290(gene="protein", state="I"),
         hsa5290(gene="protein", state="A"), [kf146, kr146, kc146])
catalyze(hsa1439(gene="protein", state="A"),
         hsa5291(gene="protein", state="I"),
         hsa5291(gene="protein", state="A"), [kf147, kr147, kc147])
catalyze(hsa1439(gene="protein", state="A"),
         hsa4893(gene="protein", state="I"),
         hsa4893(gene="protein", state="A"), [kf148, kr148, kc148])
catalyze(hsa1439(gene="protein", state="A"),
         hsa5293(gene="protein", state="I"),
         hsa5293(gene="protein", state="A"), [kf149, kr149, kc149])
catalyze(hsa1439(gene="protein", state="A"),
         hsa5294(gene="protein", state="I"),
         hsa5294(gene="protein", state="A"), [kf150, kr150, kc150])
catalyze(hsa1439(gene="protein", state="A"),
         hsa5295(gene="protein", state="I"),
         hsa5295(gene="protein", state="A"), [kf151, kr151, kc151])
catalyze(hsa1439(gene="protein", state="A"),
         hsa5296(gene="protein", state="I"),
         hsa5296(gene="protein", state="A"), [kf152, kr152, kc152])
catalyze(hsa1439(gene="protein", state="A"),
         hsa3845(gene="protein", state="I"),
         hsa3845(gene="protein", state="A"), [kf153, kr153, kc153])
catalyze(hsa1649(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf154, kr154, kc154])
catalyze(hsa1649(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf155, kr155, kc155])
catalyze(hsa1508(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf156, kr156, kc156])
catalyze(hsa1508(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf157, kr157, kc157])
catalyze(hsa1508(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf158, kr158, kc158])
catalyze(hsa1508(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf159, kr159, kc159])
catalyze(hsa1508(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf160, kr160, kc160])
catalyze(hsa1508(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf161, kr161, kc161])
catalyze(hsa1508(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf162, kr162, kc162])
catalyze(cpdC05981(bf=None), hsa8503(gene="protein", state="I"),
         hsa8503(gene="protein", state="A"), [kf163, kr163, kc163])
catalyze(hsa598(gene="protein", state="A"), hsa581(gene="protein", state="A"),
         hsa581(gene="protein", state="I"), [kf164, kr164, kc164])
catalyze(hsa598(gene="protein", state="A"), hsa578(gene="protein", state="A"),
         hsa578(gene="protein", state="I"), [kf165, kr165, kc165])
catalyze(hsa596(gene="protein", state="A"), hsa581(gene="protein", state="A"),
         hsa581(gene="protein", state="I"), [kf166, kr166, kc166])
catalyze(hsa596(gene="protein", state="A"), hsa578(gene="protein", state="A"),
         hsa578(gene="protein", state="I"), [kf167, kr167, kc167])
catalyze(hsa3562(gene="protein", state="A"),
         hsa1439(gene="protein", state="I"),
         hsa1439(gene="protein", state="A"), [kf168, kr168, kc168])
catalyze(hsa3562(gene="protein", state="A"),
         hsa3563(gene="protein", state="I"),
         hsa3563(gene="protein", state="A"), [kf169, kr169, kc169])
catalyze(hsa3563(gene="protein", state="A"),
         hsa23533(gene="protein", state="I"),
         hsa23533(gene="protein", state="A"), [kf170, kr170, kc170])
catalyze(hsa3563(gene="protein", state="A"),
         hsa8503(gene="protein", state="I"),
         hsa8503(gene="protein", state="A"), [kf171, kr171, kc171])
catalyze(hsa3563(gene="protein", state="A"),
         hsa3265(gene="protein", state="I"),
         hsa3265(gene="protein", state="A"), [kf172, kr172, kc172])
catalyze(hsa3563(gene="protein", state="A"),
         hsa5290(gene="protein", state="I"),
         hsa5290(gene="protein", state="A"), [kf173, kr173, kc173])
catalyze(hsa3563(gene="protein", state="A"),
         hsa5291(gene="protein", state="I"),
         hsa5291(gene="protein", state="A"), [kf174, kr174, kc174])
catalyze(hsa3563(gene="protein", state="A"),
         hsa4893(gene="protein", state="I"),
         hsa4893(gene="protein", state="A"), [kf175, kr175, kc175])
catalyze(hsa3563(gene="protein", state="A"),
         hsa5293(gene="protein", state="I"),
         hsa5293(gene="protein", state="A"), [kf176, kr176, kc176])
catalyze(hsa3563(gene="protein", state="A"),
         hsa5294(gene="protein", state="I"),
         hsa5294(gene="protein", state="A"), [kf177, kr177, kc177])
catalyze(hsa3563(gene="protein", state="A"),
         hsa5295(gene="protein", state="I"),
         hsa5295(gene="protein", state="A"), [kf178, kr178, kc178])
catalyze(hsa3563(gene="protein", state="A"),
         hsa5296(gene="protein", state="I"),
         hsa5296(gene="protein", state="A"), [kf179, kr179, kc179])
catalyze(hsa3563(gene="protein", state="A"),
         hsa3845(gene="protein", state="I"),
         hsa3845(gene="protein", state="A"), [kf180, kr180, kc180])
catalyze(hsa841(gene="protein", state="A"), hsa836(gene="protein", state="I"),
         hsa836(gene="protein", state="A"), [kf181, kr181, kc181])
catalyze(hsa841(gene="protein", state="A"), hsa839(gene="protein", state="I"),
         hsa839(gene="protein", state="A"), [kf182, kr182, kc182])
catalyze(hsa841(gene="protein", state="A"), hsa840(gene="protein", state="I"),
         hsa840(gene="protein", state="A"), [kf183, kr183, kc183])
catalyze(hsa840(gene="protein", state="A"), hsa142(gene="protein", state="A"),
         hsa142(gene="protein", state="I"), [kf184, kr184, kc184])
catalyze(hsa840(gene="protein", state="A"), hsa143(gene="protein", state="A"),
         hsa143(gene="protein", state="I"), [kf185, kr185, kc185])
catalyze(hsa840(gene="protein", state="A"),
         hsa10039(gene="protein", state="A"),
         hsa10039(gene="protein", state="I"), [kf186, kr186, kc186])
catalyze(hsa840(gene="protein", state="A"),
         hsa10038(gene="protein", state="A"),
         hsa10038(gene="protein", state="I"), [kf187, kr187, kc187])
catalyze(hsa840(gene="protein", state="A"), hsa71(gene="protein", state="A"),
         hsa71(gene="protein", state="I"), [kf188, kr188, kc188])
catalyze(hsa840(gene="protein", state="A"), hsa60(gene="protein", state="A"),
         hsa60(gene="protein", state="I"), [kf189, kr189, kc189])
catalyze(hsa842(gene="protein", state="A"), hsa836(gene="protein", state="I"),
         hsa836(gene="protein", state="A"), [kf190, kr190, kc190])
catalyze(hsa842(gene="protein", state="A"), hsa839(gene="protein", state="I"),
         hsa839(gene="protein", state="A"), [kf191, kr191, kc191])
catalyze(hsa842(gene="protein", state="A"), hsa142(gene="protein", state="A"),
         hsa142(gene="protein", state="I"), [kf192, kr192, kc192])
catalyze(hsa842(gene="protein", state="A"), hsa143(gene="protein", state="A"),
         hsa143(gene="protein", state="I"), [kf193, kr193, kc193])
catalyze(hsa842(gene="protein", state="A"),
         hsa10039(gene="protein", state="A"),
         hsa10039(gene="protein", state="I"), [kf194, kr194, kc194])
catalyze(hsa842(gene="protein", state="A"),
         hsa10038(gene="protein", state="A"),
         hsa10038(gene="protein", state="I"), [kf195, kr195, kc195])
catalyze(hsa842(gene="protein", state="A"), hsa840(gene="protein", state="I"),
         hsa840(gene="protein", state="A"), [kf196, kr196, kc196])
catalyze(hsa1520(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf197, kr197, kc197])
catalyze(hsa1520(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf198, kr198, kc198])
catalyze(hsa1520(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf199, kr199, kc199])
catalyze(hsa1520(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf200, kr200, kc200])
catalyze(hsa1520(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf201, kr201, kc201])
catalyze(hsa1520(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf202, kr202, kc202])
catalyze(hsa1520(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf203, kr203, kc203])
catalyze(hsa4893(gene="protein", state="A"),
         hsa5894(gene="protein", state="I"),
         hsa5894(gene="protein", state="A"), [kf204, kr204, kc204])
catalyze(hsa1075(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf205, kr205, kc205])
catalyze(hsa1075(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf206, kr206, kc206])
catalyze(hsa1075(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf207, kr207, kc207])
catalyze(hsa1075(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf208, kr208, kc208])
catalyze(hsa1075(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf209, kr209, kc209])
catalyze(hsa1075(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf210, kr210, kc210])
catalyze(hsa1075(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf211, kr211, kc211])
catalyze(cpdC05981(bf=None), hsa5170(gene="protein", state="I"),
         hsa5170(gene="protein", state="A"), [kf212, kr212, kc212])
catalyze(hsa472(gene="protein", state="A"), hsa7157(gene="protein", state="I"),
         hsa7157(gene="protein", state="A"), [kf213, kr213, kc213])
catalyze(hsa5170(gene="protein", state="A"), hsa208(gene="protein", state="I"),
         hsa208(gene="protein", state="A"), [kf214, kr214, kc214])
catalyze(hsa5170(gene="protein", state="A"), hsa207(gene="protein", state="I"),
         hsa207(gene="protein", state="A"), [kf215, kr215, kc215])
catalyze(hsa5170(gene="protein", state="A"),
         hsa10000(gene="protein", state="I"),
         hsa10000(gene="protein", state="A"), [kf216, kr216, kc216])
catalyze(hsa3002(gene="protein", state="A"),
         hsa7277(gene="protein", state="A"),
         hsa7277(gene="protein", state="I"), [kf217, kr217, kc217])
catalyze(hsa3002(gene="protein", state="A"), hsa841(gene="protein", state="I"),
         hsa841(gene="protein", state="A"), [kf218, kr218, kc218])
catalyze(hsa3002(gene="protein", state="A"), hsa840(gene="protein", state="I"),
         hsa840(gene="protein", state="A"), [kf219, kr219, kc219])
catalyze(hsa3002(gene="protein", state="A"),
         hsa1676(gene="protein", state="A"),
         hsa1676(gene="protein", state="I"), [kf220, kr220, kc220])
catalyze(hsa3002(gene="protein", state="A"),
         hsa4000(gene="protein", state="A"),
         hsa4000(gene="protein", state="I"), [kf221, kr221, kc221])
catalyze(hsa3002(gene="protein", state="A"),
         hsa4001(gene="protein", state="A"),
         hsa4001(gene="protein", state="I"), [kf222, kr222, kc222])
catalyze(hsa3002(gene="protein", state="A"),
         hsa10039(gene="protein", state="A"),
         hsa10039(gene="protein", state="I"), [kf223, kr223, kc223])
catalyze(hsa3002(gene="protein", state="A"),
         hsa10038(gene="protein", state="A"),
         hsa10038(gene="protein", state="I"), [kf224, kr224, kc224])
catalyze(hsa3002(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf225, kr225, kc225])
catalyze(hsa3002(gene="protein", state="A"),
         hsa79861(gene="protein", state="A"),
         hsa79861(gene="protein", state="I"), [kf226, kr226, kc226])
catalyze(hsa3002(gene="protein", state="A"),
         hsa84823(gene="protein", state="A"),
         hsa84823(gene="protein", state="I"), [kf227, kr227, kc227])
catalyze(hsa3002(gene="protein", state="A"), hsa143(gene="protein", state="A"),
         hsa143(gene="protein", state="I"), [kf228, kr228, kc228])
catalyze(hsa3002(gene="protein", state="A"), hsa142(gene="protein", state="A"),
         hsa142(gene="protein", state="I"), [kf229, kr229, kc229])
catalyze(hsa3002(gene="protein", state="A"),
         hsa51807(gene="protein", state="A"),
         hsa51807(gene="protein", state="I"), [kf230, kr230, kc230])
catalyze(hsa3002(gene="protein", state="A"),
         hsa7846(gene="protein", state="A"),
         hsa7846(gene="protein", state="I"), [kf231, kr231, kc231])
catalyze(hsa3002(gene="protein", state="A"), hsa836(gene="protein", state="I"),
         hsa836(gene="protein", state="A"), [kf232, kr232, kc232])
catalyze(hsa3002(gene="protein", state="A"),
         hsa7278(gene="protein", state="A"),
         hsa7278(gene="protein", state="I"), [kf233, kr233, kc233])
catalyze(hsa3002(gene="protein", state="A"),
         hsa4170(gene="protein", state="A"),
         hsa4170(gene="protein", state="I"), [kf234, kr234, kc234])
catalyze(hsa3002(gene="protein", state="A"),
         hsa84790(gene="protein", state="A"),
         hsa84790(gene="protein", state="I"), [kf235, kr235, kc235])
catalyze(hsa3002(gene="protein", state="A"),
         hsa113457(gene="protein", state="A"),
         hsa113457(gene="protein", state="I"), [kf236, kr236, kc236])
catalyze(hsa3002(gene="protein", state="A"),
         hsa112714(gene="protein", state="A"),
         hsa112714(gene="protein", state="I"), [kf237, kr237, kc237])
catalyze(hsa3002(gene="protein", state="A"),
         hsa10376(gene="protein", state="A"),
         hsa10376(gene="protein", state="I"), [kf238, kr238, kc238])
catalyze(hsa5595(gene="protein"), hsa596(gene="off"),
         hsa596(gene="on", state="A"), [kf239, kr239, kc239])
Rule("hsa596_expression_239",
     hsa596(bf=None, gene="on", state="A") >> hsa596(bf=None, gene="on",
                                                     state="A") + hsa596(
         bf=None, gene="protein", state="A"), kr239)
catalyze(hsa5594(gene="protein"), hsa596(gene="off"),
         hsa596(gene="on", state="A"), [kf240, kr240, kc240])
Rule("hsa596_expression_240",
     hsa596(bf=None, gene="on", state="A") >> hsa596(bf=None, gene="on",
                                                     state="A") + hsa596(
         bf=None, gene="protein", state="A"), kr240)
catalyze(cpdC05981(bf=None), hsa5290(gene="protein", state="I"),
         hsa5290(gene="protein", state="A"), [kf241, kr241, kc241])
catalyze(cpdC05981(bf=None), hsa5291(gene="protein", state="I"),
         hsa5291(gene="protein", state="A"), [kf242, kr242, kc242])
catalyze(cpdC05981(bf=None), hsa5293(gene="protein", state="I"),
         hsa5293(gene="protein", state="A"), [kf243, kr243, kc243])
catalyze(cpdC05981(bf=None), hsa5294(gene="protein", state="I"),
         hsa5294(gene="protein", state="A"), [kf244, kr244, kc244])
catalyze(cpdC05981(bf=None), hsa5296(gene="protein", state="I"),
         hsa5296(gene="protein", state="A"), [kf245, kr245, kc245])
catalyze(hsa4170(gene="protein", state="A"),
         hsa10018(gene="protein", state="A"),
         hsa10018(gene="protein", state="I"), [kf246, kr246, kc246])
catalyze(hsa10018(gene="protein", state="A"),
         hsa596(gene="protein", state="A"), hsa596(gene="protein", state="I"),
         [kf247, kr247, kc247])
catalyze(hsa10018(gene="protein", state="A"),
         hsa581(gene="protein", state="I"), hsa581(gene="protein", state="A"),
         [kf248, kr248, kc248])
catalyze(hsa10018(gene="protein", state="A"),
         hsa598(gene="protein", state="A"), hsa598(gene="protein", state="I"),
         [kf249, kr249, kc249])
catalyze(hsa10000(gene="protein", state="A"),
         hsa3551(gene="protein", state="I"),
         hsa3551(gene="protein", state="A"), [kf250, kr250, kc250])
catalyze(hsa10000(gene="protein", state="A"),
         hsa1147(gene="protein", state="I"),
         hsa1147(gene="protein", state="A"), [kf251, kr251, kc251])
catalyze(hsa10000(gene="protein", state="A"),
         hsa8517(gene="protein", state="I"),
         hsa8517(gene="protein", state="A"), [kf252, kr252, kc252])
catalyze(hsa10000(gene="protein", state="A"),
         hsa572(gene="protein", state="A"), hsa572(gene="protein", state="I"),
         [kf253, kr253, kc253])
catalyze(hsa56616(gene="protein", state="A"),
         hsa329(gene="protein", state="A"), hsa329(gene="protein", state="I"),
         [kf254, kr254, kc254])
catalyze(hsa56616(gene="protein", state="A"),
         hsa331(gene="protein", state="A"), hsa331(gene="protein", state="I"),
         [kf255, kr255, kc255])
catalyze(hsa56616(gene="protein", state="A"),
         hsa330(gene="protein", state="A"), hsa330(gene="protein", state="I"),
         [kf256, kr256, kc256])
catalyze(hsa56616(gene="protein", state="A"),
         hsa332(gene="protein", state="A"), hsa332(gene="protein", state="I"),
         [kf257, kr257, kc257])
catalyze(hsa332(gene="protein", state="A"), hsa836(gene="protein", state="A"),
         hsa836(gene="protein", state="I"), [kf258, kr258, kc258])
catalyze(hsa332(gene="protein", state="A"), hsa840(gene="protein", state="A"),
         hsa840(gene="protein", state="I"), [kf259, kr259, kc259])
catalyze(hsa332(gene="protein", state="A"), hsa842(gene="protein", state="A"),
         hsa842(gene="protein", state="I"), [kf260, kr260, kc260])
catalyze(hsa331(gene="protein", state="A"), hsa836(gene="protein", state="A"),
         hsa836(gene="protein", state="I"), [kf261, kr261, kc261])
catalyze(hsa331(gene="protein", state="A"), hsa840(gene="protein", state="A"),
         hsa840(gene="protein", state="I"), [kf262, kr262, kc262])
catalyze(hsa331(gene="protein", state="A"), hsa842(gene="protein", state="A"),
         hsa842(gene="protein", state="I"), [kf263, kr263, kc263])
catalyze(hsa330(gene="protein", state="A"), hsa836(gene="protein", state="A"),
         hsa836(gene="protein", state="I"), [kf264, kr264, kc264])
catalyze(hsa330(gene="protein", state="A"), hsa840(gene="protein", state="A"),
         hsa840(gene="protein", state="I"), [kf265, kr265, kc265])
catalyze(hsa330(gene="protein", state="A"), hsa842(gene="protein", state="A"),
         hsa842(gene="protein", state="I"), [kf266, kr266, kc266])
catalyze(hsa824(gene="protein", state="A"),
         hsa100506742(gene="protein", state="I"),
         hsa100506742(gene="protein", state="A"), [kf267, kr267, kc267])
catalyze(hsa823(gene="protein", state="A"),
         hsa100506742(gene="protein", state="I"),
         hsa100506742(gene="protein", state="A"), [kf268, kr268, kc268])
catalyze(hsa8837(gene="protein", state="A"), hsa841(gene="protein", state="A"),
         hsa841(gene="protein", state="I"), [kf269, kr269, kc269])
catalyze(hsa8837(gene="protein", state="A"), hsa843(gene="protein", state="A"),
         hsa843(gene="protein", state="I"), [kf270, kr270, kc270])
catalyze(hsa1616(gene="protein", state="A"),
         hsa4217(gene="protein", state="I"),
         hsa4217(gene="protein", state="A"), [kf271, kr271, kc271])
catalyze(hsa5605(gene="protein", state="A"),
         hsa5595(gene="protein", state="I"),
         hsa5595(gene="protein", state="A"), [kf272, kr272, kc272])
catalyze(hsa5605(gene="protein", state="A"),
         hsa5594(gene="protein", state="I"),
         hsa5594(gene="protein", state="A"), [kf273, kr273, kc273])
catalyze(hsa5604(gene="protein", state="A"),
         hsa5595(gene="protein", state="I"),
         hsa5595(gene="protein", state="A"), [kf274, kr274, kc274])
catalyze(hsa5604(gene="protein", state="A"),
         hsa5594(gene="protein", state="I"),
         hsa5594(gene="protein", state="A"), [kf275, kr275, kc275])
catalyze(hsa5602(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf276, kr276, kc276])
catalyze(hsa5602(gene="protein", state="A"), hsa572(gene="protein", state="I"),
         hsa572(gene="protein", state="A"), [kf277, kr277, kc277])
catalyze(hsa5602(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf278, kr278, kc278])
catalyze(hsa5602(gene="protein", state="A"),
         hsa3725(gene="protein", state="I"),
         hsa3725(gene="protein", state="A"), [kf279, kr279, kc279])
catalyze(hsa5602(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf280, kr280, kc280])
catalyze(hsa5602(gene="protein", state="A"),
         hsa2353(gene="protein", state="I"),
         hsa2353(gene="protein", state="A"), [kf281, kr281, kc281])
catalyze(hsa5601(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf282, kr282, kc282])
catalyze(hsa5601(gene="protein", state="A"), hsa572(gene="protein", state="I"),
         hsa572(gene="protein", state="A"), [kf283, kr283, kc283])
catalyze(hsa5601(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf284, kr284, kc284])
catalyze(hsa5601(gene="protein", state="A"),
         hsa3725(gene="protein", state="I"),
         hsa3725(gene="protein", state="A"), [kf285, kr285, kc285])
catalyze(hsa5601(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf286, kr286, kc286])
catalyze(hsa5601(gene="protein", state="A"),
         hsa2353(gene="protein", state="I"),
         hsa2353(gene="protein", state="A"), [kf287, kr287, kc287])
catalyze(hsa8743(gene="protein", state="A"),
         hsa8793(gene="protein", state="I"),
         hsa8793(gene="protein", state="A"), [kf288, kr288, kc288])
catalyze(hsa8743(gene="protein", state="A"),
         hsa8797(gene="protein", state="I"),
         hsa8797(gene="protein", state="A"), [kf289, kr289, kc289])
catalyze(hsa8743(gene="protein", state="A"),
         hsa8794(gene="protein", state="I"),
         hsa8794(gene="protein", state="A"), [kf290, kr290, kc290])
catalyze(hsa8743(gene="protein", state="A"),
         hsa8795(gene="protein", state="I"),
         hsa8795(gene="protein", state="A"), [kf291, kr291, kc291])
catalyze(hsa4914(gene="protein", state="A"),
         hsa23533(gene="protein", state="I"),
         hsa23533(gene="protein", state="A"), [kf292, kr292, kc292])
catalyze(hsa4914(gene="protein", state="A"),
         hsa8503(gene="protein", state="I"),
         hsa8503(gene="protein", state="A"), [kf293, kr293, kc293])
catalyze(hsa4914(gene="protein", state="A"),
         hsa3845(gene="protein", state="I"),
         hsa3845(gene="protein", state="A"), [kf294, kr294, kc294])
catalyze(hsa4914(gene="protein", state="A"),
         hsa5290(gene="protein", state="I"),
         hsa5290(gene="protein", state="A"), [kf295, kr295, kc295])
catalyze(hsa4914(gene="protein", state="A"),
         hsa5291(gene="protein", state="I"),
         hsa5291(gene="protein", state="A"), [kf296, kr296, kc296])
catalyze(hsa4914(gene="protein", state="A"),
         hsa4893(gene="protein", state="I"),
         hsa4893(gene="protein", state="A"), [kf297, kr297, kc297])
catalyze(hsa4914(gene="protein", state="A"),
         hsa5293(gene="protein", state="I"),
         hsa5293(gene="protein", state="A"), [kf298, kr298, kc298])
catalyze(hsa4914(gene="protein", state="A"),
         hsa5294(gene="protein", state="I"),
         hsa5294(gene="protein", state="A"), [kf299, kr299, kc299])
catalyze(hsa4914(gene="protein", state="A"),
         hsa5295(gene="protein", state="I"),
         hsa5295(gene="protein", state="A"), [kf300, kr300, kc300])
catalyze(hsa4914(gene="protein", state="A"),
         hsa5296(gene="protein", state="I"),
         hsa5296(gene="protein", state="A"), [kf301, kr301, kc301])
catalyze(hsa4914(gene="protein", state="A"),
         hsa3265(gene="protein", state="I"),
         hsa3265(gene="protein", state="A"), [kf302, kr302, kc302])
catalyze(hsa9451(gene="protein", state="A"),
         hsa1965(gene="protein", state="I"),
         hsa1965(gene="protein", state="A"), [kf303, kr303, kc303])
catalyze(cpdC05981(bf=None), hsa5295(gene="protein", state="I"),
         hsa5295(gene="protein", state="A"), [kf304, kr304, kc304])
catalyze(hsa329(gene="protein", state="A"), hsa836(gene="protein", state="A"),
         hsa836(gene="protein", state="I"), [kf305, kr305, kc305])
catalyze(hsa329(gene="protein", state="A"), hsa840(gene="protein", state="A"),
         hsa840(gene="protein", state="I"), [kf306, kr306, kc306])
catalyze(hsa329(gene="protein", state="A"), hsa842(gene="protein", state="A"),
         hsa842(gene="protein", state="I"), [kf307, kr307, kc307])
catalyze(hsa5599(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf308, kr308, kc308])
catalyze(hsa5599(gene="protein", state="A"), hsa572(gene="protein", state="I"),
         hsa572(gene="protein", state="A"), [kf309, kr309, kc309])
catalyze(hsa5599(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf310, kr310, kc310])
catalyze(hsa5599(gene="protein", state="A"),
         hsa3725(gene="protein", state="I"),
         hsa3725(gene="protein", state="A"), [kf311, kr311, kc311])
catalyze(hsa5599(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf312, kr312, kc312])
catalyze(hsa5599(gene="protein", state="A"),
         hsa2353(gene="protein", state="I"),
         hsa2353(gene="protein", state="A"), [kf313, kr313, kc313])
catalyze(hsa836(gene="protein", state="A"), hsa71(gene="protein", state="A"),
         hsa71(gene="protein", state="I"), [kf314, kr314, kc314])
catalyze(hsa836(gene="protein", state="A"), hsa142(gene="protein", state="A"),
         hsa142(gene="protein", state="I"), [kf315, kr315, kc315])
catalyze(hsa836(gene="protein", state="A"), hsa143(gene="protein", state="A"),
         hsa143(gene="protein", state="I"), [kf316, kr316, kc316])
catalyze(hsa836(gene="protein", state="A"),
         hsa10039(gene="protein", state="A"),
         hsa10039(gene="protein", state="I"), [kf317, kr317, kc317])
catalyze(hsa836(gene="protein", state="A"),
         hsa10038(gene="protein", state="A"),
         hsa10038(gene="protein", state="I"), [kf318, kr318, kc318])
catalyze(hsa836(gene="protein", state="A"), hsa6709(gene="protein", state="A"),
         hsa6709(gene="protein", state="I"), [kf319, kr319, kc319])
catalyze(hsa836(gene="protein", state="A"), hsa6708(gene="protein", state="A"),
         hsa6708(gene="protein", state="I"), [kf320, kr320, kc320])
catalyze(hsa836(gene="protein", state="A"), hsa1676(gene="protein", state="A"),
         hsa1676(gene="protein", state="I"), [kf321, kr321, kc321])
catalyze(hsa836(gene="protein", state="A"), hsa60(gene="protein", state="A"),
         hsa60(gene="protein", state="I"), [kf322, kr322, kc322])
catalyze(hsa5414(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf323, kr323, kc323])
catalyze(hsa5414(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf324, kr324, kc324])
catalyze(hsa5414(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf325, kr325, kc325])
catalyze(hsa5414(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf326, kr326, kc326])
catalyze(hsa1521(gene="protein", state="A"), hsa637(gene="protein", state="I"),
         hsa637(gene="protein", state="A"), [kf327, kr327, kc327])
catalyze(hsa1521(gene="protein", state="A"), hsa598(gene="protein", state="A"),
         hsa598(gene="protein", state="I"), [kf328, kr328, kc328])
catalyze(hsa1521(gene="protein", state="A"), hsa596(gene="protein", state="A"),
         hsa596(gene="protein", state="I"), [kf329, kr329, kc329])
catalyze(hsa1521(gene="protein", state="A"), hsa329(gene="protein", state="A"),
         hsa329(gene="protein", state="I"), [kf330, kr330, kc330])
catalyze(hsa1521(gene="protein", state="A"), hsa332(gene="protein", state="A"),
         hsa332(gene="protein", state="I"), [kf331, kr331, kc331])
catalyze(hsa1521(gene="protein", state="A"), hsa331(gene="protein", state="A"),
         hsa331(gene="protein", state="I"), [kf332, kr332, kc332])
catalyze(hsa1521(gene="protein", state="A"), hsa330(gene="protein", state="A"),
         hsa330(gene="protein", state="I"), [kf333, kr333, kc333])
catalyze(hsa839(gene="protein", state="A"), hsa60(gene="protein", state="A"),
         hsa60(gene="protein", state="I"), [kf334, kr334, kc334])
catalyze(hsa839(gene="protein", state="A"), hsa71(gene="protein", state="A"),
         hsa71(gene="protein", state="I"), [kf335, kr335, kc335])
catalyze(hsa839(gene="protein", state="A"),
         hsa84823(gene="protein", state="A"),
         hsa84823(gene="protein", state="I"), [kf336, kr336, kc336])
catalyze(hsa839(gene="protein", state="A"), hsa4000(gene="protein", state="A"),
         hsa4000(gene="protein", state="I"), [kf337, kr337, kc337])
catalyze(hsa839(gene="protein", state="A"), hsa4001(gene="protein", state="A"),
         hsa4001(gene="protein", state="I"), [kf338, kr338, kc338])
catalyze(hsa7157(gene="protein"), hsa637(gene="off"),
         hsa637(gene="on", state="A"), [kf339, kr339, kc339])
Rule("hsa637_expression_339",
     hsa637(bf=None, gene="on", state="A") >> hsa637(bf=None, gene="on",
                                                     state="A") + hsa637(
         bf=None, gene="protein", state="A"), kr339)
catalyze(hsa7157(gene="protein"), hsa317(gene="off"),
         hsa317(gene="on", state="A"), [kf340, kr340, kc340])
Rule("hsa317_expression_340",
     hsa317(bf=None, gene="on", state="A") >> hsa317(bf=None, gene="on",
                                                     state="A") + hsa317(
         bf=None, gene="protein", state="A"), kr340)
catalyze(hsa7157(gene="protein"), hsa8793(gene="off"),
         hsa8793(gene="on", state="A"), [kf341, kr341, kc341])
Rule("hsa8793_expression_341",
     hsa8793(bf=None, gene="on", state="A") >> hsa8793(bf=None, gene="on",
                                                       state="A") + hsa8793(
         bf=None, gene="protein", state="A"), kr341)
catalyze(hsa7157(gene="protein"), hsa578(gene="off"),
         hsa578(gene="on", state="A"), [kf342, kr342, kc342])
Rule("hsa578_expression_342",
     hsa578(bf=None, gene="on", state="A") >> hsa578(bf=None, gene="on",
                                                     state="A") + hsa578(
         bf=None, gene="protein", state="A"), kr342)
catalyze(hsa7157(gene="protein"), hsa8794(gene="off"),
         hsa8794(gene="on", state="A"), [kf343, kr343, kc343])
Rule("hsa8794_expression_343",
     hsa8794(bf=None, gene="on", state="A") >> hsa8794(bf=None, gene="on",
                                                       state="A") + hsa8794(
         bf=None, gene="protein", state="A"), kr343)
catalyze(hsa7157(gene="protein"), hsa581(gene="off"),
         hsa581(gene="on", state="A"), [kf344, kr344, kc344])
Rule("hsa581_expression_344",
     hsa581(bf=None, gene="on", state="A") >> hsa581(bf=None, gene="on",
                                                     state="A") + hsa581(
         bf=None, gene="protein", state="A"), kr344)
catalyze(hsa7157(gene="protein"), hsa27113(gene="off"),
         hsa27113(gene="on", state="A"), [kf345, kr345, kc345])
Rule("hsa27113_expression_345",
     hsa27113(bf=None, gene="on", state="A") >> hsa27113(bf=None, gene="on",
                                                         state="A") + hsa27113(
         bf=None, gene="protein", state="A"), kr345)
catalyze(hsa7157(gene="protein"), hsa63970(gene="off"),
         hsa63970(gene="on", state="A"), [kf346, kr346, kc346])
Rule("hsa63970_expression_346",
     hsa63970(bf=None, gene="on", state="A") >> hsa63970(bf=None, gene="on",
                                                         state="A") + hsa63970(
         bf=None, gene="protein", state="A"), kr346)
catalyze(hsa7157(gene="protein"), hsa5366(gene="off"),
         hsa5366(gene="on", state="A"), [kf347, kr347, kc347])
Rule("hsa5366_expression_347",
     hsa5366(bf=None, gene="on", state="A") >> hsa5366(bf=None, gene="on",
                                                       state="A") + hsa5366(
         bf=None, gene="protein", state="A"), kr347)
catalyze(hsa7157(gene="protein"), hsa355(gene="off"),
         hsa355(gene="on", state="A"), [kf348, kr348, kc348])
Rule("hsa355_expression_348",
     hsa355(bf=None, gene="on", state="A") >> hsa355(bf=None, gene="on",
                                                     state="A") + hsa355(
         bf=None, gene="protein", state="A"), kr348)
catalyze(hsa7157(gene="protein"), hsa8795(gene="off"),
         hsa8795(gene="on", state="A"), [kf349, kr349, kc349])
Rule("hsa8795_expression_349",
     hsa8795(bf=None, gene="on", state="A") >> hsa8795(bf=None, gene="on",
                                                       state="A") + hsa8795(
         bf=None, gene="protein", state="A"), kr349)
catalyze(hsa7157(gene="protein"), hsa8797(gene="off"),
         hsa8797(gene="on", state="A"), [kf350, kr350, kc350])
Rule("hsa8797_expression_350",
     hsa8797(bf=None, gene="on", state="A") >> hsa8797(bf=None, gene="on",
                                                       state="A") + hsa8797(
         bf=None, gene="protein", state="A"), kr350)
catalyze(hsa7157(gene="protein"), hsa55367(gene="off"),
         hsa55367(gene="on", state="A"), [kf351, kr351, kc351])
Rule("hsa55367_expression_351",
     hsa55367(bf=None, gene="on", state="A") >> hsa55367(bf=None, gene="on",
                                                         state="A") + hsa55367(
         bf=None, gene="protein", state="A"), kr351)
catalyze(hsa4217(gene="protein", state="A"),
         hsa5601(gene="protein", state="I"),
         hsa5601(gene="protein", state="A"), [kf352, kr352, kc352])
catalyze(hsa4217(gene="protein", state="A"),
         hsa5599(gene="protein", state="I"),
         hsa5599(gene="protein", state="A"), [kf353, kr353, kc353])
catalyze(hsa4217(gene="protein", state="A"),
         hsa5602(gene="protein", state="I"),
         hsa5602(gene="protein", state="A"), [kf354, kr354, kc354])
catalyze(hsa1965(gene="protein", state="A"), hsa468(gene="protein", state="I"),
         hsa468(gene="protein", state="A"), [kf355, kr355, kc355])
catalyze(hsa3265(gene="protein", state="A"),
         hsa5894(gene="protein", state="I"),
         hsa5894(gene="protein", state="A"), [kf356, kr356, kc356])
