#!/usr/bin/env python

# Run TargetScan's context score script (targetscan_60_context_scores.pl) on the output of TargetScan's base
# script (targetscan_60.pl). The script takes as input the id and sequence of an siRNA, as well as the directory
# of seeds scanned with TargetScan's base script and a translation file from tarnscripts to genes, and outputs
# TargetScan context score predictions to an output folder. Remember, we cannot use aligned UTRs, since siRNAs are
# not endogenous entities.
#
# Author: Fabian Schmich (fabian.schmich@bsse.ethz.ch)
#

# Libraries
import sys
import os.path
import optparse
import tempfile
from numpy import mean

# Paths
CS_SCRIPT = "/path/to/targetscan_60_context_scores.pl"
UTR_FILE = "./data/hg19_ucsc_3p.txt"
SCANNED_SEEDS_DIR = "/path/to/seeds/"
TRANSLATION = "./data/refseq.txt"

class Target:
    "Class to represent a siRNA target gene"
    def __init__(self, id, name=None, scores=[], transcripts=[]):
        self.id = id
        self.name = name
        self.scores = scores
        self.transcripts = transcripts
    def __hash__(self):
        return hash(self.id)
    def __eq__(self, other):
        return self.id == other.id
    def __repr__(self):
        if self.name is not None: return "\t".join([str(self.id), self.name])
        else: return str(self.id)
    def strextd(self):
        return "\t".join([str(self.id), self.name, ";".join(self.transcripts), ";".join([str(e) for e in self.scores]), str(mean(self.scores))]) 

# Takes an siRNAs and produces the input files required by context_score script,
# "miRNA file" and "PredictedTargets" file by matching the sRNA seed to pre-scanned
# seeds.
def prepare(sirna_id, sirna_seq, species="9606"):
    # Mature sequence file
    mature_sequence_file = tempfile.NamedTemporaryFile("w", prefix="seq_", delete=False)
    mature_sequence_file.write("\t".join(["miRNA_family_ID", "Species_ID", "MiRBase_ID", "Mature_sequence"])+"\n")
    mature_sequence_file.write("\t".join([sirna_id, str(species), sirna_id, sirna_seq])+"\n")
    mature_sequence_file.close()
    # Seed predictions file
    ts_predictions_file = tempfile.NamedTemporaryFile("w", prefix="targets_", delete=False)
    seedsf = open(SCANNED_SEEDS_DIR + "/" + sirna_seq[1:8] + ".txt", "r")
    ts_predictions_file.write(seedsf.readline()) # header
    for line in seedsf:
        lsp = line.split("\t")
        lsp[1] = sirna_id
        ts_predictions_file.write("\t".join(lsp))
    seedsf.close()
    ts_predictions_file.close()
    return(mature_sequence_file.name, ts_predictions_file.name)


def predict(mat_seq_file_name, ts_pred_file_name):
    contextplus_score_file = tempfile.NamedTemporaryFile("w", prefix="cps_", delete=True)
    contextplus_score_file.close()
    olddir = os.getcwd()
    os.chdir(os.path.dirname(CS_SCRIPT))
    from subprocess import call
    call([CS_SCRIPT, mat_seq_file_name, UTR_FILE, ts_pred_file_name, contextplus_score_file.name])
    os.chdir(olddir)
    return (contextplus_score_file.name)


def parse(csfile):    
    # Hash translation
    pos_tr = 0
    pos_gid = 1
    pos_gn = 2
    translf = open(TRANSLATION, "r")
    tr = dict()
    for line in translf:
        lsp = line.split("\t")
        lsp = [el.strip() for el in lsp]
        if lsp[pos_tr] not in tr:
            tr[lsp[pos_tr]] = [[lsp[pos_gid], lsp[pos_gn]]]
        else:
            tr[lsp[pos_tr]].append([lsp[pos_gid], lsp[pos_gn]])
    translf.close()
    # Parse
    pos_trname = 0
    pos_score = 11
    transcripts = dict()
    csf = open(csfile, "r")
    csf.readline() # header
    for line in csf:
        lsp = line.split("\t")
        if lsp[pos_score].strip() == "too_close": 
            continue
        else:
            score = round(float(lsp[pos_score]), 3)
        trname = lsp[pos_trname].strip()
        if trname in transcripts:
            transcripts[trname].append(score)
        else:
            transcripts[trname] = [score]
    csf.close()
    # Translate
    targets = []
    for trkey in transcripts.keys():
        trname = trkey
        if "." in trname: 
            trname = trname[:trname.index(".")] # Strip everything after "."
        if trname in tr:
            for [gid, gn] in tr[trname]:
                curgene = Target(id=gid, name=gn)
                if curgene in targets:
                    g_ind = targets.index(curgene)
                    targets[g_ind].scores.append(sum(transcripts[trkey]))
                    targets[g_ind].transcripts.append(trname)
                else:
                    curgene.scores = [sum(transcripts[trkey])]
                    curgene.transcripts = [trname]
                    targets.append(curgene)
        else:
            print("WARN: cannot translate transcript %s" % trname)
    return(targets)


def main():
    # Parse Arguments
    usage = "Usage: %prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option(
        "--id",
        action="store",
        type="string",
        dest="id",
        help="siRNA identifier"
    )
    parser.add_option(
        "--seq",
        action="store",
        type="string",
        dest="seq",
        help="siRNA antisense sequence w/o overhangs"
    )
    parser.add_option(
        "--out",
        action="store",
        type="string",
        dest="out",
        help="Output file path"
    )
    (options, args) = parser.parse_args(sys.argv)
    # Parse input argumeeans
    if options.id is None or options.seq is None or options.out is None:
        print "ERROR: Missing siRNA identifier or sequence"
        sys.exit(0)
    print "%s: %s" % (options.id, options.seq)
    # Prepare input files
    [mat_seq_file, ts_pred_file] = prepare(options.id, options.seq)
    # Predict context plus scores
    contextplus_score_file = predict(mat_seq_file, ts_pred_file)
    # Parse context plus scores and translate transcripts to genes
    targets = parse(contextplus_score_file)
    # Write out targets
    outf = open(options.out, "w")
    outf.write("GeneID\tGeneName\tTranscripts\tCPS\tCPSmean\n")
    for target in targets:
        outf.write(target.strextd()+"\n")    
    outf.close()
    # Remove temp files
    for f in [mat_seq_file, ts_pred_file, contextplus_score_file]: 
        os.remove(f)

if __name__ == "__main__":
    main()
