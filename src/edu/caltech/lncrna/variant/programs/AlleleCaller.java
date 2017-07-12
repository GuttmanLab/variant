package edu.caltech.lncrna.variant.programs;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.caltech.lncrna.bio.alignment.Alignment;
import edu.caltech.lncrna.bio.alignment.ChromosomeName;
import edu.caltech.lncrna.bio.alignment.CoordinateSpace;
import edu.caltech.lncrna.bio.alignment.SingleReadAlignment;
import edu.caltech.lncrna.bio.io.BamWriter;
import edu.caltech.lncrna.bio.io.SingleReadBamParser;
import edu.caltech.lncrna.bio.sequence.Base;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public final class AlleleCaller {

    private Path inputVcfPath;
    private Path inputBamPath;
    private Path outputBamPath1;
    private Path outputBamPath2;
    private Path outputBamPathAmbg;
    private Path outputBamPathConf;
    private String genotype1;
    private String genotype2;
    private boolean diffChromNotation;
    
    private long totalCount;
    private long gt1Count;
    private long gt2Count;
    private long ambgCount;
    private long confCount;
    
    private static final String BAM = ".bam";
    private static final int PROGRESS_DENOMINATOR = 100000;
    
    private static final String TRANSLATION_ERROR = "A problem occured when " +
            "translating between different chromosome notations. This " +
            "program only supports \"chr1\"-to-\"1\" conversions, and vice-" +
            "versa. (An exception is made for \"chrM\"-to-\"MT\".)";
    
    private static final Logger LOGGER =
            Logger.getLogger(AlleleCaller.class.getName());
    
    public static void main(String[] args) throws IOException {
        (new AlleleCaller())
            .parseArgs(args)
            .run()
            .logOutput();
    }
    
    private AlleleCaller() {
        LOGGER.setLevel(Level.INFO);
    }
    
    private AlleleCaller parseArgs(String[] args) {
        Options options = new Options();

        options.addOption(Option.builder()
                .longOpt("vcf")
                .desc("input VCF file")
                .argName("VCF")
                .hasArg()
                .required()
                .build());

        options.addOption(Option.builder()
                .longOpt("bam")
                .desc("input BAM file of reads to decouple")
                .argName("FILE")
                .hasArg()
                .required()
                .build());
        
        options.addOption(Option.builder()
                .longOpt("gt1")
                .desc("name of genotype 1")
                .argName("NAME")
                .hasArg()
                .required()
                .build());
        
        options.addOption(Option.builder()
                .longOpt("gt2")
                .desc("name of genotype 1")
                .argName("NAME")
                .hasArg()
                .required()
                .build());
        
        options.addOption(Option.builder()
                .longOpt("convert")
                .desc("name of chromosomes in VCF file and BAM file do not " +
                        "agree (e.g., \"chr1\" and \"1\")")
                .hasArg(false)
                .required(false)
                .build());
        
        options.addOption(Option.builder()
                .longOpt("debug")
                .desc("set logging level to debug")
                .required(false)
                .build());

        try {
            CommandLineParser parser = new DefaultParser();
            CommandLine cmd = parser.parse(options, args);
            inputVcfPath = Paths.get(cmd.getOptionValue("vcf"));            
            String bamString = cmd.getOptionValue("bam");
            inputBamPath = Paths.get(bamString);
            genotype1 = cmd.getOptionValue("gt1");
            genotype2 = cmd.getOptionValue("gt2");
            diffChromNotation = cmd.hasOption("convert");
            if (cmd.hasOption("debug")) {
                LOGGER.setLevel(Level.FINEST);
            }
        } catch (ParseException e) {
            printHelp(options);
            System.exit(-1);
        }
        
        String inputBamString = inputBamPath.toString();
        
        outputBamPath1 = Paths.get(inputBamString.substring(0,
                inputBamString.length() - BAM.length()) + "." + genotype1 +
                BAM);
        outputBamPath2 = Paths.get(inputBamString.substring(0,
                inputBamString.length() - BAM.length()) + "." + genotype2 +
                BAM);
        outputBamPathAmbg = Paths.get(inputBamString.substring(0,
                inputBamString.length() - BAM.length()) + ".ambiguous" +
                BAM);
        outputBamPathConf = Paths.get(inputBamString.substring(0,
                inputBamString.length() - BAM.length()) + ".conflicting" +
                BAM);
        
        checkVcfFile();

        return this;
    }
    
    private void checkVcfFile() {
        try (VCFFileReader vcfReader = new VCFFileReader(inputVcfPath.toFile())) {
            List<String> gts = vcfReader.getFileHeader().getGenotypeSamples();
            if (!gts.contains(genotype1)) {
                LOGGER.log(Level.SEVERE,
                        "VCF file does not contain genotype " + genotype1);
                System.exit(-1);
            }
            if (!gts.contains(genotype2)) {
                LOGGER.log(Level.SEVERE,
                        "VCF file does not contain genotype " + genotype2);
                System.exit(-1);
            }
        } catch (TribbleException | RuntimeIOException ex) {
            LOGGER.log(Level.SEVERE, "Error opening VCF file.", ex);
            System.exit(-1);
        }
    }
    
    /**
     * Creates and prints a help menu to the console.
     * @param opts - the options to include in the help menu
     */
    private void printHelp(Options opts) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.setDescPadding(0);
        String header = "\n";
        String footer = "\n";
        formatter.printHelp("java -jar VariantCaller.jar", header,
                opts, footer, true);
    }
    
    private AlleleCaller run() {

        CoordinateSpace header = new CoordinateSpace(inputBamPath);
        
        try (SingleReadBamParser bp = new SingleReadBamParser(inputBamPath);
             VCFFileReader vcfReader = new VCFFileReader(inputVcfPath.toFile());
             BamWriter out1 = new BamWriter(outputBamPath1, header);
             BamWriter out2 = new BamWriter(outputBamPath2, header);
             BamWriter outAmbg = new BamWriter(outputBamPathAmbg, header);
             BamWriter outConf = new BamWriter(outputBamPathConf, header)) {
            Iterator<SingleReadAlignment> reads = bp.getAlignmentIterator();
            while (reads.hasNext()) {
                SingleReadAlignment read = reads.next();
                Classification c = classify(read, vcfReader);
                
                LOGGER.finest(read.getName() + " classified as " + c);
                
                switch (c) {
                case VAR1:
                    out1.writeSamRecord(read);
                    gt1Count++;
                    break;
                case VAR2:
                    out2.writeSamRecord(read);
                    gt2Count++;
                    break;
                case CONF:
                    outConf.writeSamRecord(read);
                    confCount++;
                    break;
                default:
                    ambgCount++;
                    outAmbg.writeSamRecord(read);
                }
                if (++totalCount % PROGRESS_DENOMINATOR == 0) {
                    LOGGER.info("Processing read " + totalCount + ".");
                }
            }
        }

        return this;
    }
    
    /**
     * Logs the final count totals.
     */
    private void logOutput() {
        LOGGER.info(gt1Count + " reads assigned to " + genotype1 + ".");
        LOGGER.info(gt2Count + " reads assigned to " + genotype2 + ".");
        LOGGER.info(ambgCount + " reads were ambiguous.");
        LOGGER.info(confCount + " reads were in conflict.");
    }
    
    /**
     * Assigns a read to a genotype by considering all of its SNPs.
     */
    private Classification classify(Alignment read, VCFFileReader vcfReader) {
        Classification rtrn = Classification.AMBG;
        String chr = diffChromNotation
                ? translate(read.getReferenceName())
                : read.getReferenceName();
        int start = read.getStart();
        int end = read.getEnd();
        CloseableIterator<VariantContext> snps = vcfReader.query(chr, start, end);
                        
        while (snps.hasNext()) {
            VariantContext snp = snps.next();
            rtrn = rtrn.and(classify(read, snp));
            
            if (rtrn.equals(Classification.CONF)) {
                return rtrn;
            }
        }

        return rtrn;
    }
    
    /**
     * Assigns a SNP to a genotype.
     */
    private Classification classify(Alignment read, VariantContext snp) {
        Genotype gt1 = snp.getGenotype(genotype1);
        Genotype gt2 = snp.getGenotype(genotype2);

        if (gt1.sameGenotype(gt2)) {
            return Classification.AMBG;
        }
        
        String gtString1 = gt1.getGenotypeString();
        String gtString2 = gt2.getGenotypeString();
        
        // Not sure how to deal with genotype strings longer than 3
        if (gtString1.length() != 3 || gtString2.length() != 3) {
            return Classification.AMBG;
        }
        
        if (gtString1.charAt(0) != gtString1.charAt(2)) {
            return Classification.AMBG;
        }
        
        if (gtString2.charAt(0) != gtString2.charAt(2)) {
            return Classification.AMBG;
        }
        
        Base base1 = Base.of(gtString1.charAt(0));
        Base base2 = Base.of(gtString2.charAt(0));
        
        // This case should have been caught when comparing genotypes
        assert !base1.equals(base2);

        Base b = read.getReadBaseFromReferencePosition(snp.getStart());
        
        if (b.equals(base1)) {
            return Classification.VAR1;
        } else if (b.equals(base2)) {
            return Classification.VAR2;
        } else {
            return Classification.CONF;
        }
    }
    
    /**
     * Translates between chromosome notations.
     * <p>
     * Converts "chr1" to "1" and vice versa. The only exception is the pair
     * "chrM" and "MT", which map to one another.
     */
    private String translate(String chr) {
        String rtrn = null;
        try {
            rtrn = ChromosomeName.MAPPING.get(chr);
            if (chr == null) {
                throw new NoSuchElementException("Chromosome " + chr +
                        " not found in the lookup table.");
            }
        } catch (NoSuchElementException ex) {
            LOGGER.log(Level.SEVERE, TRANSLATION_ERROR, ex);
            System.exit(-1);
        }
        return rtrn;
    }
    
    /**
     * An enumeration of how to classify a read given a SNP
     * <ul>
     * <li>VAR1: belongs to genotype 1
     * <li>VAR2: belongs to genotype 2
     * <li>AMBG: no difference between either genotype
     * <li>CONF: conflicting, contains SNPs from both genotypes
     * </ul>
     */
    private enum Classification {
        VAR1() {
            public Classification and(Classification other) {
                switch (other) {
                case VAR1:
                case AMBG:
                    return this;
                case VAR2:
                case CONF:
                default:
                    return Classification.CONF;
                }
            }
        },
        
        VAR2() {
            public Classification and(Classification other) {
                switch (other) {
                case VAR2:
                case AMBG:
                    return this;
                case VAR1:
                case CONF:
                default:
                    return Classification.CONF;
                }
            }            
        },

        AMBG() {
            public Classification and(Classification other) {
                return other;
            }
        },

        CONF() {
            public Classification and(Classification other) {
                return this;
            }
        };
        
        /**
         * Combines this classification with another.
         * <p>
         * Classifications are initially determined for each SNP position in a
         * read. Some of these SNPs may be ambiguous and provide no
         * information; other SNPs may conflict with one another. This method
         * combines the classifications from two SNPs (and should be chained to
         * deal with any number of SNPs) to produce a final classification for
         * the read as a whole.
         */
        public abstract Classification and(Classification other);
    }
}