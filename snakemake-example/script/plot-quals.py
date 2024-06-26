import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

quals = [record.qual for record in VariantFile("output/calls/all.vcf")]
plt.hist(quals)

plt.savefig("output/plots/quals.svg")
