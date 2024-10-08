'''
The main purpose of this script is to create read length histogram plots for 
read sequencing data.
Input data can be given as one or multiple of:
-compressed, standard or streamed fastq file
-compressed, standard or streamed fastq file, with
 additional information added by albacore or MinKNOW
-a bam file
-a summary file generated by albacore
'''

from os import path
import logging
import nanomath
import numpy as np
from scipy import stats
import nanoplot.utils as utils
import nanoplot.report as report
from nanoget import get_input
from nanoplot.filteroptions import filter_and_transform_data
from nanoplot.version import __version__
import nanoplotter
import pickle
import sys
from nanoplot.NanoPlot import make_stats
import beepy as bp

def hist_main():
    '''
    Organization function
    -setups logging
    -gets inputdata
    -calls plotting function
    '''
    settings, args = utils.get_args()
    try:
        utils.make_output_dir(args.outdir)
        utils.init_logs(args)
        # args.format = nanoplotter.check_valid_format(args.format)
        if args.pickle:
            datadf = pickle.load(open(args.pickle, 'rb'))
        elif args.feather:
            from nanoget import combine_dfs
            from pandas import read_feather
            datadf = combine_dfs([read_feather(p) for p in args.feather], method="simple")
        else:
            sources = {
                "fastq": args.fastq,
                "bam": args.bam,
                "cram": args.cram,
                "fastq_rich": args.fastq_rich,
                "fastq_minimal": args.fastq_minimal,
                "summary": args.summary,
                "fasta": args.fasta,
                "ubam": args.ubam,
            }
            datadf = get_input(
                source=[n for n, s in sources.items() if s][0],
                files=[f for f in sources.values() if f][0],
                threads=args.threads,
                readtype=args.readtype,
                combine="simple",
                barcoded=args.barcoded,
                huge=args.huge,
                keep_supp=not (args.no_supplementary))
        if args.store:
            pickle.dump(
                obj=datadf,
                file=open(settings["path"] + "NanoPlot-data.pickle", 'wb'))
        if args.raw:
            datadf.to_csv(settings["path"] + "NanoPlot-data.tsv.gz",
                          sep="\t",
                          index=False,
                          compression="gzip")

        settings["statsfile"] = [make_stats(datadf, settings, suffix="", tsv_stats=args.tsv_stats)]
        datadf, settings = filter_and_transform_data(datadf, settings)
        # if settings["filtered"]:  # Bool set when filter was applied in filter_and_transform_data()
        #     settings["statsfile"].append(
        #         make_stats(datadf[datadf["length_filter"]], settings,
        #                    suffix="_post_filtering", tsv_stats=args.tsv_stats)
        #     )

        # if args.barcoded:
        #     main_path = settings["path"]
        #     barcodes = list(datadf["barcode"].unique())
        #     plots = []
        #     for barc in barcodes:
        #         logging.info("Processing {}".format(barc))
        #         dfbarc = datadf[datadf["barcode"] == barc]
        #         if len(dfbarc) > 5:
        #             settings["title"] = barc
        #             settings["path"] = path.join(args.outdir, args.prefix + barc + "_")
        #             plots.append(report.BarcodeTitle(barc))
        #             plots.extend(
        #                 make_plots(dfbarc, settings)
        #             )
        #         else:
        #             sys.stderr.write("Found barcode {} less than 5x, ignoring...\n".format(barc))
        #             logging.info("Found barcode {} less than 5 times, ignoring".format(barc))
        #     settings["path"] = main_path
        # else:
        #     plots = make_plots(datadf, settings)
        # make_report(plots, settings)
        plots = make_hist_plots(datadf, settings, args)

        logging.info("Finished!")
        bp.beep(sound=1)

    except Exception as e:
        logging.error(e, exc_info=True)
        print("\n\n\nIf you read this then NanoPlot {} has crashed :-(".format(__version__))
        print("Please try updating NanoPlot and see if that helps...\n")
        print("If not, please report this issue at https://github.com/wdecoster/NanoPlot/issues")
        print("If you could include the log file that would be really helpful.")
        print("Thanks!\n\n\n")
        raise

def make_hist_plots(datadf, settings, args):
    '''
    Call plotting functions from nanoplotter
    settings["lengths_pointer"] is a column in the DataFrame specifying which lengths to use
    '''
    color = nanoplotter.check_valid_color(settings["color"])
    colormap = nanoplotter.check_valid_colormap(settings["colormap"])

    plotdict = {type: settings["plots"].count(type) for type in ["kde", "hex", "dot", 'pauvre']}
    if "hex" in settings["plots"]:
        print(
            "WARNING: hex as part of --plots has been deprecated and will be ignored. To get the hex output, rerun with --legacy hex.")

    if settings["legacy"] is None:
        plotdict_legacy = {}
    else:
        plotdict_legacy = {plot: settings["legacy"].count(plot) for plot in ["kde", "hex", "dot"]}
    plots = []

    subdf = utils.subsample_datasets(datadf)
    if settings["N50"]:
        n50 = nanomath.get_N50(np.sort(datadf["lengths"]))
    else:
        n50 = None

    # Make histogram plots with outliers
    plots.extend(
        nanoplotter.length_plots(
            array=datadf[datadf["length_filter"]]["lengths"].astype('uint64'),
            name="Read length",
            path=settings["path"],
            n50=n50,
            color=color,
            title=settings["title"],
            figformat=settings["format"])
    )
    print("Created length plots with outliers")
    logging.info("Created length plots with outliers")

    # Make histogram plots without outliers
    settings['drop_outliers'] = True
    datadf, settings = filter_and_transform_data(datadf, settings)
    settings["statsfile"].append(
        make_stats(datadf[datadf["length_filter"]], settings,
                   suffix="_post_filtering", tsv_stats=args.tsv_stats)
    )

    plots.extend(
        nanoplotter.length_plots(
            array=datadf[datadf["length_filter"]]["lengths"].astype('uint64'),
            name="Read length without outliers",
            path=settings["path"],
            n50=n50,
            color=color,
            title=settings["title"],
            figformat=settings["format"])
    )
    print("Created length plots without outliers")
    logging.info("Created length plots without outliers")
    return plots

if __name__ == "__main__":
    hist_main()