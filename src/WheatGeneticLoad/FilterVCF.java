/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import daxing.applets.ScriptMethods;
import format.position.ChrPos;
import format.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.IOFileFormat;
import utils.IOUtils;
import utils.PArrayUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class FilterVCF {

    public FilterVCF() {
        //this.statVcfCoverage();
        //this.subsetCovevsSDvsPV();
        //this.statVcfPValue();
        //this.subsetCovevsSDvsPV();
        //this.gethighdensitySNPPos();
        //this.addGroup();

        //new ScriptMethods().getCellDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/abd/000_chr1A-7A.ABDgenome.depth.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrAsub.ABDgenome");
        //new ScriptMethods().getCellDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/000_chr1B-7B.ABDgenome.depth.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrBsub.ABDgenome/");
        //new ScriptMethods().getCellDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/000_chr1D-7D.ABDgenome.depth.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrDsub.ABDgenome/");
        //this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrAsub.ABDgenome/position","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.65.depthVSsd.txt.gz", 140);
        //this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrAsub.ABDgenome/position","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.75.depthVSsd.txt.gz", 226);
        //this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrAsub.ABDgenome/position","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.85.depthVSsd.txt.gz", 348);
        //this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrAsub.ABDgenome/position","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.80.depthVSsd.txt.gz", 282);
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1A-7A.ABDgenome.depth_5208sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrAsub.ABDgenome.depthVSsd.addGroup.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1A-7A.ABDgenome.depth_5208sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrAsub.ABDgenome.depthVSsd.addGroup.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1A-7A.ABDgenome.depth_5208sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrAsub.ABDgenome.depthVSsd.addGroup.bin100_0.85.txt");
        //this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.80.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1A-7A.ABDgenome.depth_5208sites.txt", "80", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrAsub.ABDgenome.depthVSsd.addGroup.bin100_0.80.txt");
        //this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrBsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.75.depthVSsd.txt.gz", 248);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrBsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.65.depthVSsd.txt.gz", 156);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrBsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.85.depthVSsd.txt.gz", 381);
//this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrBsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.80.depthVSsd.txt.gz", 309);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrDsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.75.depthVSsd.txt.gz", 129);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrDsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.65.depthVSsd.txt.gz", 75);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrDsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.85.depthVSsd.txt.gz", 233);
        //this.randomTxt();
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1B-7B.ABDgenome.depth_7000sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrBsub.ABDgenome.depthVSsd.addGroup.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1B-7B.ABDgenome.depth_7000sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrBsub.ABDgenome.depthVSsd.addGroup.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1B-7B.ABDgenome.depth_7000sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrBsub.ABDgenome.depthVSsd.addGroup.bin100_0.85.txt");       
        //this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.80.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1B-7B.ABDgenome.depth_4996sites.txt", "80", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrBsub.ABDgenome.depthVSsd.addGroup.bin100_0.80.txt");       
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1D-7D.ABDgenome.depth_7000sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrDsub.ABDgenome.depthVSsd.addGroup.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1D-7D.ABDgenome.depth_7000sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrDsub.ABDgenome.depthVSsd.addGroup.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1D-7D.ABDgenome.depth_7000sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/abd/chrDsub.ABDgenome.depthVSsd.addGroup.bin100_0.85.txt");
//        new ScriptMethods().getCellDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/001_chr1A-7A.ABgenome.depth.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/ab/chrAsub.ABgenome");
//        new ScriptMethods().getCellDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/001_chr1B-7B.ABgenome.depth.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/ab/chrBsub.ABgenome");
//        new ScriptMethods().getCellDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/002_chr1D-7D.Dgenome.depth.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/d");
//        
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/ab/chrAsub.ABgenome/position", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrAsub/chrAsub_AB_bin100_0.65.depthVSsd.txt.gz", 61);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/ab/chrAsub.ABgenome/position", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrAsub/chrAsub_AB_bin100_0.75.depthVSsd.txt.gz", 102);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/ab/chrAsub.ABgenome/position", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrAsub/chrAsub_AB_bin100_0.85.depthVSsd.txt.gz", 173);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/ab/chrBsub.ABgenome/position", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrBsub/chrBsub_AB_bin100_0.65.depthVSsd.txt.gz", 80);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/ab/chrBsub.ABgenome/position", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrBsub/chrBsub_AB_bin100_0.75.depthVSsd.txt.gz", 129);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/ab/chrBsub.ABgenome/position", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrBsub/chrBsub_AB_bin100_0.85.depthVSsd.txt.gz", 206);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/d/position", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/d/chrDgenome_bin100_0.65.depthVSsd.txt.gz", 203);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/d/position", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/d/chrDgenome_bin100_0.75.depthVSsd.txt.gz", 232);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/d/position", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/d/chrDgenome_bin100_0.85.depthVSsd.txt.gz", 562);
        //       this.randomTxt();
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrAsub/chrAsub_AB_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1A-7A.ABgenome.depth_4955sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/ab/chrAsub.ABgenome.depthVSsd.addGroup.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrAsub/chrAsub_AB_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1A-7A.ABgenome.depth_4955sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/ab/chrAsub.ABgenome.depthVSsd.addGroup.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrAsub/chrAsub_AB_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1A-7A.ABgenome.depth_4955sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/ab/chrAsub.ABgenome.depthVSsd.addGroup.bin100_0.85.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrBsub/chrBsub_AB_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1B-7B.ABgenome.depth_4908sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/ab/chrBsub.ABgenome.depthVSsd.addGroup.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrBsub/chrBsub_AB_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1B-7B.ABgenome.depth_4908sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/ab/chrBsub.ABgenome.depthVSsd.addGroup.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrBsub/chrBsub_AB_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1B-7B.ABgenome.depth_4908sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/ab/chrBsub.ABgenome.depthVSsd.addGroup.bin100_0.85.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/d/chrDgenome_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/d/002_chr1D-7D.Dgenome.depth_5049sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/d/chrDgenome.depthVSsd.addGroup.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/d/chrDgenome_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/d/002_chr1D-7D.Dgenome.depth_5049sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/d/chrDgenome.depthVSsd.addGroup.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/d/chrDgenome_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/d/002_chr1D-7D.Dgenome.depth_5049sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addGroup/d/chrDgenome.depthVSsd.addGroup.bin100_0.85.txt");
////
        //this.mergePosList();
        
        
    }

    /**
     * 将chr001.ABDgenome 和 chr001.ABgenome 合并起来，生成一个合并后的文件
     */
    public void mergePosList(String inFileS1, String inFileS2, String outfileS ) {
//        String inFileS1 = "/Users/Aoyue/Documents/chr001.ABDgenome.10000lines.vcf.gz";
//        String inFileS2 = "/Users/Aoyue/Documents/chr001.ABgenome.10000lines.vcf.gz";
//        String outfileS = "/Users/Aoyue/Documents/chr001_PosAllele.txt";
        String f1 = new File(inFileS1).getName();
        String f2 = new File(inFileS2).getName();
        String[] alleles = {"A", "C", "G", "T", "D", "I"};
        Arrays.sort(alleles);
        double[] fre = new double[6];
        try {
            int chr1 = Integer.MIN_VALUE;
            int chr2 = Integer.MIN_VALUE;
            int taxaNum1 = Integer.MIN_VALUE;
            int taxaNum2 = Integer.MIN_VALUE;
            TIntArrayList posList1 = new TIntArrayList();
            List<String> referList1 = new ArrayList<>();
            List<String> altList1 = new ArrayList<>();
            List<String> altDepthList1 = new ArrayList<>();
            TIntArrayList posList2 = new TIntArrayList();
            List<String> referList2 = new ArrayList<>();
            List<String> altList2 = new ArrayList<>();
            List<String> altDepthList2 = new ArrayList<>();
            BufferedReader br = IOUtils.getTextGzipReader(inFileS1);
            String temp = null;
            while ((temp = br.readLine()).startsWith("##")) {
            };
            taxaNum1 = temp.split("\t").length - 9;
            String[] tem = null;
            int cnt1 = 0;
            while ((temp = br.readLine()) != null) {
                cnt1++;
                temp = temp.substring(0, 50);
                tem = temp.split("\t");
                chr1 = Integer.parseInt(tem[0]);
                posList1.add(Integer.parseInt(tem[1]));
                referList1.add(tem[3]);
                altList1.add(tem[4]);
                altDepthList1.add(tem[7].split(";")[1].replace("AD=", ""));
            }
            br.close();
            System.out.println(f1 + " contains " + cnt1 + " variants.");
            br = IOUtils.getTextGzipReader(inFileS2);
            temp = null;
            while ((temp = br.readLine()).startsWith("##")) {
            };
            taxaNum2 = temp.split("\t").length - 9;
            double weight1 = (double) taxaNum1 / (taxaNum1 + taxaNum2);
            double weight2 = (double) taxaNum2 / (taxaNum1 + taxaNum2);
            tem = null;
            int cnt2 = 0;
            while ((temp = br.readLine()) != null) {
                cnt2++;
                temp = temp.substring(0, 150);
                tem = temp.split("\t");
                chr2 = Integer.parseInt(tem[0]);
                posList2.add(Integer.parseInt(tem[1]));
                referList2.add(tem[3]);
                altList2.add(tem[4]);
                altDepthList2.add(tem[7].split(";")[1].replace("AD=", ""));
            }
            System.out.println(f2 + " contains " + cnt2 + " variants.");

            if (chr1 != chr2) {
                System.out.println("Wrong input files! Program quits.");
                System.exit(0);
            }
            //************************ 总共的variants数目，F1独有的数目，F2独有的数目，F1 F2共有的数目

            int totalvariants = 0;
            int intf1 = 0;
            int intf2 = 0;
            int shared = 0;
            System.out.println("totalVariants\t" + f1 + "\t" + f2 + "\tSharedVariants");

            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tRef\tAlt");
            bw.newLine();
            TIntHashSet mergedPosSet = new TIntHashSet(posList1);
            mergedPosSet.addAll(posList2);
            int[] mergedPos = mergedPosSet.toArray();
            Arrays.sort(mergedPos);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < mergedPos.length; i++) {
                totalvariants++;
                sb = new StringBuilder();
                int index1 = posList1.binarySearch(mergedPos[i]);
                int index2 = posList2.binarySearch(mergedPos[i]);
                if (index1 < 0 && index2 > -1) {
                    if (altList2.get(index2).length() > 3) {
                        continue;
                    }
                    intf2++;
                    sb.append(chr1).append("\t").append(posList2.get(index2)).append("\t").append(referList2.get(index2)).append("\t").append(altList2.get(index2));
                } else if (index1 > -1 && index2 < 0) {
                    if (altList1.get(index1).length() > 3) {
                        continue;
                    }
                    intf1++;
                    sb.append(chr1).append("\t").append(posList1.get(index1)).append("\t").append(referList1.get(index1)).append("\t").append(altList1.get(index1));
                } else {
                    shared++;
                    for (int j = 0; j < fre.length; j++) { //初始化fre,使之都等于-1
                        fre[j] = -1;
                    }
                    //先处理ABD群体
                    tem = altList1.get(index1).split(",");// 总共含有的alt数目 ，在I的循环里，搜索 posList对应的index,根据index找到alt的信息，放入tem中
                    String[] fretem = altDepthList1.get(index1).split(","); //AD的深度数组
                    double[] depth = new double[fretem.length]; //
                    double[] fre1 = new double[depth.length];
                    double sum = 0;
                    for (int j = 0; j < depth.length; j++) { //根据AD的深度个数，求总的depth
                        depth[j] = Integer.parseInt(fretem[j]);
                        sum += depth[j];
                    }
                    for (int j = 0; j < depth.length; j++) { //求各个Alt的等位基因频率
                        fre1[j] = depth[j] / sum;
                    }
                    for (int j = 0; j < tem.length; j++) { //fre指的是"A", "C", "G", "T", "D", "I"
                        int index = Arrays.binarySearch(alleles, tem[j]);//第0，1个AD，在 ACGTDI中的index搜索
                        fre[index] = fre1[j+1] * weight1; //复制fre1[j+1]到 fre中去， 权重1等于 在ABD中群体的个数 除以 在ABD中群体的个数加上AB群体的个数之和
                    }
                    //再处理AB群体
                    tem = altList2.get(index2).split(",");
                    fretem = altDepthList2.get(index2).split(",");
                    depth = new double[fretem.length];
                    double[] fre2 = new double[depth.length];
                    sum = 0;
                    for (int j = 0; j < depth.length; j++) {
                        depth[j] = Integer.parseInt(fretem[j]);
                        sum += depth[j];
                    }
                    for (int j = 0; j < depth.length; j++) {
                        fre2[j] = depth[j] / sum;
                    }
                    for (int j = 0; j < tem.length; j++) {
                        int index = Arrays.binarySearch(alleles, tem[j]);
                        if (fre[index] < 0) {
                            fre[index] = fre2[j+1] * weight2; //如果没搜到，说明在ABD群体中没有发现，只在群体2中发现；
                        } else {
                            fre[index] = fre[index] + fre2[j+1] * weight2; //说明在ABD和AB群体中都有，2者相加。
                        }
                    }

                    int[] indices = PArrayUtils.getIndexByDescendingValue(fre);
                    sb.append(chr1).append("\t").append(posList1.get(index1)).append("\t").append(referList1.get(index1)).append("\t");
                    for (int j = 0; j < 2; j++) {
                        if (fre[indices[j]] > 0) {
                            sb.append(alleles[indices[j]]).append(",");
                        } else {
                            break;
                        }
                    }
                    sb.deleteCharAt(sb.length() - 1);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            
            System.out.println(mergedPos.length + "\t" + (intf1 + shared) + "\t" + (intf2 + shared) + "\t" + shared);
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * 这里采用单线程 对抽样的vcf文件进行每个位点深度和sd进行统计和Pvalue获取，制成一个表格，chr pos averageDepth SD
     * PValue 注意表格不要以#开头，否则被注释，看不到
     *
     */
    public void statVcfDepth_SD_PValue_singlethread(String infileS, String outfileS) {

        BufferedReader br = null;
        if (infileS.endsWith(".vcf")) {
            br = IOUtils.getTextReader(infileS);
        } else if (infileS.endsWith(".vcf.gz")) {
            br = IOUtils.getTextGzipReader(infileS);
        }
        BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
        String temp = null;
        int cnt = 0;
        try {
            while ((temp = br.readLine()).startsWith("##")) {
            }
            List<String> linetaxa = PStringUtils.fastSplit(temp, "\t");
            bw.write(linetaxa.get(0).replaceFirst("#", "") + "\t" + linetaxa.get(1) + "\t" + "AverageDepth\tSD\tPValue");
            bw.newLine();
            String[] taxa = new String[linetaxa.size() - 9];
            cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++; // 对snp开始计数
                if (cnt % 1000000 == 0) {
                    System.out.println(String.valueOf(cnt) + " lines");
                }
                TDoubleArrayList depthList = new TDoubleArrayList();
                TDoubleArrayList PValueList = new TDoubleArrayList();
                List<String> l = new ArrayList<>();
                l = PStringUtils.fastSplit(temp, "\t");
                //排除错误的那一行
                String chr = l.get(0);
                String pos = l.get(1);
//                if ((chr.equals("2")) && (pos.equals("94099978"))) {
//                    continue;
//                }
//                if ((chr.equals("34")) && (pos.equals("263599921"))) {
//                    continue;
//                }
                if ((chr.equals("33")) && (pos.equals("391299785"))) {
                    continue;
                }
                String pvalue = l.get(7).split("PV=")[1].split(";")[0];
                for (int i = 0; i < taxa.length; i++) {
                    String genoS = l.get(i + 9);
                    if (genoS.startsWith(".")) {
                        depthList.add(0);
                        continue;
                    }
                    List<String> ll = PStringUtils.fastSplit(genoS, ":");
                    List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                    int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                    depthList.add(depth);
                }
                double[] dep = depthList.toArray();
                DescriptiveStatistics d = new DescriptiveStatistics(dep);
                double relativeMean = d.getMean();
                double sd = d.getStandardDeviation();
                //计算完毕，接下来开始写入文件
                StringBuilder sb = new StringBuilder();
                sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(String.format("%.6f", relativeMean)).append("\t").append(String.format("%.6f", sd))
                        .append("\t").append(pvalue);
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println(temp);
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println(infileS + " is calculated well done");
    }

    //1.从DepthDB中随机抽取5000个位点后，根据上一步生成的0.65 0.75 0.85 的总pos信息，为随机抽样的位点添加group， 0.65分组为6 0.75分组为7； 0.85分组为8
    public void addGroupforsubsetDepthDB(String cumuFileS, String depthfileS, String groupnum, String outfileS) {
        try {
            //1.先建立分组的库
            List<ChrPos> l = new ArrayList();
            BufferedReader br = IOUtils.getTextGzipReader(cumuFileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                String chr = PStringUtils.fastSplit(temp).get(0);
                String pos = PStringUtils.fastSplit(temp).get(1);
                l.add(new ChrPos(Short.valueOf(chr), Integer.valueOf(pos)));
            }
            br.close();
            Collections.sort(l);

            //2.读取要添加分组的文件，进行chrpos判断，为其加上分组；
            br = IOUtils.getTextReader(depthfileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine() + "\tGroup"); //读取表头
            bw.newLine();
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                String chr = PStringUtils.fastSplit(temp).get(0);
                String pos = PStringUtils.fastSplit(temp).get(1);
                int index = Collections.binarySearch(l, new ChrPos(Short.valueOf(chr), Integer.valueOf(pos)));
                if (index < 0) {
                    sb.append(temp).append("\t0");
                } else {
                    sb.append(temp).append("\t").append(groupnum);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    //随机抽样txt文本
    public void randomTxt() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/ab";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/d";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/d/";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/abd/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd";

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);

/////////////////////////////////////////////////////////////////////////////
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS;
                BufferedReader br = null;
                infileS = f.getAbsolutePath();
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                //开始数行数
                int totallines = 0;
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    totallines++;
                }
                br.close();
                double ratio = (double) 5000 / totallines; //注意一定要在5000千加上 强制类型转换，不然不能得出小数

                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String outfileS = new File(outfileDirS, f.getName().replaceFirst("depth.txt.gz", "depth_7000sites.txt")).getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write(br.readLine());
                bw.newLine();
                while ((temp = br.readLine()) != null) {
                    double r = Math.random();
                    if (r > ratio) {
                        continue;
                    }
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
/////////////////////////////////////////////////////////////////////////////

    }

    public void getHighCumulativePos(String inputDirOfPosition, String outputFile, int num) {
        File[] files = IOUtils.listRecursiveFiles(new File(inputDirOfPosition));
        List<String> chrPosList = new ArrayList<>();
        BufferedReader[] brs = new BufferedReader[num + 1];
        try {
            for (int i = 0; i < brs.length; i++) {
                brs[i] = IOUtils.getTextReader(files[i].getAbsolutePath());
                String line;
                brs[i].readLine();
                List<String> lines = new ArrayList<>();
                while ((line = brs[i].readLine()) != null) {
                    lines.add(line);
                }
                System.out.println(i + "\t" + lines.size());
                chrPosList.addAll(lines);
                brs[i].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        BufferedWriter bw = IOUtils.getTextGzipWriter(outputFile);
//        int[] randoms= ArrayTool.getRandomNonrepetitionArray(5000, 0, chrPosList.size());
//        Arrays.sort(randoms);
        try {
            bw.write("CHR" + "\t" + "POS" + "\t" + "AverageDepth" + "\t" + "SD" + "\n");
//            int index;
            for (int i = 0; i < chrPosList.size(); i++) {
//                index=Arrays.binarySearch(randoms, i);
//                if (index<0) continue;
                bw.write(chrPosList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 先建立5个数据库， 再将5000个位点读入，最后进行5个index的判断，根据5个index判断进行分组
     */
    public void addGroup() {
        String infileDirS = "";
        String dbS = "";
        String outfileS = "";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        try {
            //①先建立5个数据库， 再将5000个位点读入，最后进行5个index的判断，根据5个index判断进行分组
            List<ChrPos>[] l = new ArrayList[fs.length];
            for (int i = 0; i < fs.length; i++) {
                l[i] = new ArrayList();

            }
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    String chr = PStringUtils.fastSplit(temp).get(0);
                    String pos = PStringUtils.fastSplit(temp).get(1);
                    l[i].add(new ChrPos(Short.valueOf(chr), Integer.valueOf(pos)));
                }
            }
            System.out.println("Database has been built");

            // ②将5000个位点读入，最后进行5个index的判断
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            bw.write(br.readLine()); //读入表头
            bw.newLine();
            br.close();

            // chrAsub.ABDgenome.depthVSsd.addGroup.bin100_0.55.txt 文件名字
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String temp = br.readLine(); //read header
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                //System.out.println(String.valueOf(fs[i].getName()) + "\t" + cnt);
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 根据库文件，输出库里的pos信息，将chr pos depth sd 作图
     */
    public void gethighdensitySNPPos() {
        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/0.55-0.95/chr1A-7A_ABD_bin100_0.95.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/000_source/abd/000_chr1A-7A.ABDgenome.depth.txt.gz";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/0.55-0.95_depthVSsd/chr1A-7A_ABD_bin100_0.95.depthVSsd.txt.gz";

        List<ChrPos> posl = new ArrayList<>();
        try {
            //先建立pos数据库
            BufferedReader br = IOUtils.getTextReader(infileS1);
            BufferedWriter bw = null;
            String temp = br.readLine();
            int cnt = 0;
            short chr;
            int pos;
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                List<String> l = PStringUtils.fastSplit(temp);
                chr = Short.parseShort(l.get(0));
                pos = Integer.parseInt(l.get(1));
                posl.add(new ChrPos(chr, pos));
                cnt++;
            }
            br.close();
            System.out.println(cnt + "  pos in the database totally");

            //找出库中的pos,并写出
            Collections.sort(posl);
            br = IOUtils.getTextGzipReader(infileS2);
            bw = IOUtils.getTextGzipWriter(outfileS1);
            bw.write(br.readLine());
            bw.newLine();
            int share = 0;
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                sb.append(temp);
                List<String> l = PStringUtils.fastSplit(temp);
                chr = Short.parseShort(l.get(0));
                pos = Integer.parseInt(l.get(1));
                int index = Collections.binarySearch(posl, new ChrPos(chr, pos));
                if (index >= 0) {
                    share++;
                    bw.write(sb.toString());
                    bw.newLine();

                }
            }

            bw.flush();
            bw.close();
            br.close();
            System.out.println(share + " pos was extracted totally ");

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 对抽样的vcf文件进行每个位点深度和sd进行统计和Pvalue获取，制成一个表格，chr pos averageDepth SD PValue
     * 注意表格不要以#开头，否则被注释，看不到
     *
     */
    public void statVcfDepth_SD_PValue(String infileDirS, String outfileDirS) {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/abd/001_subset/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/abd/";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/ab/001_subsetVCF/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/ab/";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/d/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/d/";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/mergeTaxa/";
        //String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/";
        //String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/mergeTaxa/";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);

        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            //infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/chr005.merge.vcf.gz";
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".depth.txt.gz")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> linetaxa = PStringUtils.fastSplit(temp, "\t");
                bw.write(linetaxa.get(0).replaceFirst("#", "") + "\t" + linetaxa.get(1) + "\t" + "AverageDepth\tSD\tPValue");
                bw.newLine();
                String[] taxa = new String[linetaxa.size() - 9];
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++; // 对snp开始计数
                    if (cnt % 1000000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }

                    TDoubleArrayList depthList = new TDoubleArrayList();
                    TDoubleArrayList PValueList = new TDoubleArrayList();
                    List<String> l = new ArrayList<>();
                    l = PStringUtils.fastSplit(temp, "\t");
                    //排除错误的那一行
                    String chr = l.get(0);
                    String pos = l.get(1);
                    if ((chr.equals("2")) && (pos.equals("94099978"))) {
                        continue;
                    }
                    String pvalue = l.get(7).split("PV=")[1].split(";")[0];
                    for (int i = 0; i < taxa.length; i++) {
                        String genoS = l.get(i + 9);
                        if (genoS.startsWith(".")) {
                            depthList.add(0);
                            continue;
                        }
                        List<String> ll = PStringUtils.fastSplit(genoS, ":");
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                        depthList.add(depth);
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    double relativeMean = d.getMean();
                    double sd = d.getStandardDeviation();
                    //计算完毕，接下来开始写入文件
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(String.format("%.6f", relativeMean)).append("\t").append(String.format("%.6f", sd))
                            .append("\t").append(pvalue);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(f.getName() + " is calculated well done");
        });
    }

    /**
     * 对已经生成的2万行结果进行随机抽样，抽出5000行进行接下来的随机分析
     */
    public void subsetCovevsSDvsPV() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/abd/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/005_pvCalSample/abd/";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/ab/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/005_pvCalSample/ab/";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/d/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/005_pvCalSample/d/";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/";

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            String outfileS = new File(outfileDirS, f.getName().replaceFirst("depth.txt.gz", "depth_7000sites.txt")).getAbsolutePath();
            RowTable<String> t = new RowTable<>(f.getAbsolutePath());
            boolean[] ifOut = new boolean[t.getRowNumber()];
            int totallines = t.getRowNumber();
            double ratio = (double) 7000 / totallines; //注意一定要在5000千加上 强制类型转换，不然不能得出小数
            for (int i = 0; i < t.getRowNumber(); i++) {
                double r = Math.random();
                if (r > ratio) {

                } else {
                    ifOut[i] = true;
                }
            }
            t.writeTextTable(outfileS, IOFileFormat.Text, ifOut);
        });
    }

    /**
     * 对抽样的vcf文件进行每个位点每个taxa的深度统计和Pvalue获取，制成一个表格，chr pos averageDepth SD PValue
     * IfPVzero taxa1Depth ..... 注意表格不要以#开头，否则被注释，看不到
     * 代码，即在之前depth和SD的基础上又加了2列，第一列是PValue 第二列是PValue是否是0的判断,如果是，那么值为1，如果不是那么值为0
     */
    public void statVcfPValue() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/abd/001_subset/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/abd/";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/ab/001_subsetVCF/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/ab/";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/d/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/d/";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/mergeTaxa/";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/004_pvCal/mergeTaxa/";

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);

        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            //infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/013_mergeTaxaVCF/merge/chr005.merge.vcf.gz";
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".depth.txt.gz")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> l = PStringUtils.fastSplit(temp, "\t");
                bw.write(l.get(0).replaceFirst("#", "") + "\t" + l.get(1) + "\t" + "AverageDepth\tSD\tPValue\tIfPVzero");
                String[] taxa = new String[l.size() - 9];
                for (int i = 0; i < taxa.length; i++) {
                    taxa[i] = l.get(i + 9);
                    bw.write("\t" + taxa[i]);
                }
                bw.newLine();
                //建立一个整型list类型的数组，每个元素是一个list,一共有 taxa.length个list
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    TDoubleArrayList depthList = new TDoubleArrayList();
                    TDoubleArrayList PValueList = new TDoubleArrayList();

                    cnt++; // 对snp开始计数
                    if (cnt % 1000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }
                    l = PStringUtils.fastSplit(temp, "\t");
                    String pvalue = l.get(7).split("PV=")[1].split(";")[0];
                    String ifPVzero = null;
                    if (pvalue.equals("0.0")) {
                        ifPVzero = "1";
                    } else {
                        ifPVzero = "0";
                    }
                    for (int i = 0; i < taxa.length; i++) {

                        String genoS = l.get(i + 9);
                        if (genoS.startsWith(".")) {
                            depthList.add(0);
                            continue;
                        }
                        List<String> ll = PStringUtils.fastSplit(genoS, ":");
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                        depthList.add(depth);
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    double relativeMean = d.getMean();
                    double sd = d.getStandardDeviation();
                    //计算完毕，接下来开始写入文件
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(String.format("%.3f", relativeMean)).append("\t").append(String.format("%.3f", sd))
                            .append("\t").append(pvalue).append("\t").append(ifPVzero);
                    for (int i = 0; i < taxa.length; i++) {
                        sb.append("\t").append(String.format("%.0f", depthList.get(i)));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(f.getName() + " is calculated well done");
        });
    }

    /**
     * 对已经生成的2万行结果进行随机抽样，抽出5000行进行接下来的随机分析
     */
    public void subsetCoveVSSD() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/abd/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/002_depthCalSample/abd/";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/ab/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/002_depthCalSample/ab/";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/d/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/002_depthCalSample/d/";

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".depth.txt.gz", ".depth_5000sites.txt")).getAbsolutePath();
            RowTable<String> t = new RowTable<>(f.getAbsolutePath());
            boolean[] ifOut = new boolean[t.getRowNumber()];
            int totallines = t.getRowNumber();
            double ratio = (double) 5000 / totallines; //注意一定要在5000千加上 强制类型转换，不然不能得出小数
            for (int i = 0; i < t.getRowNumber(); i++) {
                double r = Math.random();
                if (r > ratio) {

                } else {
                    ifOut[i] = true;
                }
            }
            t.writeTextTable(outfileS, IOFileFormat.Text, ifOut);
        });
    }

    /**
     * 对抽样的vcf文件进行每个位点每个taxa的深度统计，制成一个表格，chr pos averageDepth SD taxa1Depth
     * ..... 注意表格不要以#开头，否则被注释，看不到
     */
    public void statVcfCoverage() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/abd/001_subset/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/abd/";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/ab/001_subsetVCF/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/ab/";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/003_subsetVCF/d/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/001_depthCal/d/";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".depth.txt.gz")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> l = PStringUtils.fastSplit(temp, "\t");
                bw.write(l.get(0).replaceFirst("#", "") + "\t" + l.get(1) + "\t" + "AverageDepth\tSD");
                String[] taxa = new String[l.size() - 9];
                for (int i = 0; i < taxa.length; i++) {
                    taxa[i] = l.get(i + 9);
                    bw.write("\t" + taxa[i]);
                }
                bw.newLine();
                //建立一个整型list类型的数组，每个元素是一个list,一共有 taxa.length个list
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    TDoubleArrayList depthList = new TDoubleArrayList();
                    cnt++; // 对snp开始计数
                    if (cnt % 1000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }
                    l = PStringUtils.fastSplit(temp, "\t");
                    for (int i = 0; i < taxa.length; i++) {

                        String genoS = l.get(i + 9);
                        if (genoS.startsWith(".")) {
                            depthList.add(0);
                            continue;
                        } //0/0:9,0:0,3,14
                        List<String> ll = PStringUtils.fastSplit(genoS, ":");
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                        depthList.add(depth);
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    double relativeMean = d.getMean(); //所有taxa的平均深度
                    double sd = d.getStandardDeviation();
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(relativeMean).append("\t").append(sd);
                    for (int i = 0; i < taxa.length; i++) {
                        sb.append("\t").append(String.format("%.0f", depthList.get(i)));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(f.getName() + " is calculated well done");
        });
    }

}
