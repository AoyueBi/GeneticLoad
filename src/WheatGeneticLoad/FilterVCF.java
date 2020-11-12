/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.CalVCF;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.table.RowTable;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import pgl.infra.pos.ChrPos;

/**
 *
 * @author Aoyue
 */
public class FilterVCF {

    public FilterVCF() {
//        this.statVcfCoverage();
//        this.subsetCovevsSDvsPV();
        //this.statVcfPValue();
//        this.subsetCovevsSDvsPV();
        //this.gethighdensitySNPPos();
        //this.addMultipleColumn();

        //new ScriptMethods().getCellDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/abd/000_chr1A-7A.ABDgenome.depth.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrAsub.ABDgenome");
        //new ScriptMethods().getCellDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/000_chr1B-7B.ABDgenome.depth.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrBsub.ABDgenome/");
        //new ScriptMethods().getCellDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/001_depthDB/000_chr1D-7D.ABDgenome.depth.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrDsub.ABDgenome/");
        //this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrAsub.ABDgenome/position","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.65.depthVSsd.txt.gz", 140);
        //this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrAsub.ABDgenome/position","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.75.depthVSsd.txt.gz", 226);
        //this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrAsub.ABDgenome/position","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.85.depthVSsd.txt.gz", 348);
        //this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrAsub.ABDgenome/position","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.80.depthVSsd.txt.gz", 282);
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1A-7A.ABDgenome.depth_5208sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrAsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1A-7A.ABDgenome.depth_5208sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrAsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1A-7A.ABDgenome.depth_5208sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrAsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.85.txt");
        //this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrAsub/chrAsub_ABD_bin100_0.80.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1A-7A.ABDgenome.depth_5208sites.txt", "80", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrAsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.80.txt");
        //this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrBsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.75.depthVSsd.txt.gz", 248);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrBsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.65.depthVSsd.txt.gz", 156);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrBsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.85.depthVSsd.txt.gz", 381);
//this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrBsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.80.depthVSsd.txt.gz", 309);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrDsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.75.depthVSsd.txt.gz", 129);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrDsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.65.depthVSsd.txt.gz", 75);
//        this.getHighCumulativePos("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/002_cellDepthDB/abd/chrDsub.ABDgenome/position/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.85.depthVSsd.txt.gz", 233);
        //this.randomTxt();
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1B-7B.ABDgenome.depth_7000sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrBsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1B-7B.ABDgenome.depth_7000sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrBsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1B-7B.ABDgenome.depth_7000sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrBsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.85.txt");
        //this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrBsub/chrBsub_ABD_bin100_0.80.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1B-7B.ABDgenome.depth_4996sites.txt", "80", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrBsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.80.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1D-7D.ABDgenome.depth_7000sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrDsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1D-7D.ABDgenome.depth_7000sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrDsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/abd/chrDsub/chrDsub_ABD_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/abd/000_chr1D-7D.ABDgenome.depth_7000sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/abd/chrDsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.85.txt");
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
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrAsub/chrAsub_AB_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1A-7A.ABgenome.depth_4955sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/ab/chrAsub.ABgenome.depthVSsd.addMultipleColumn.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrAsub/chrAsub_AB_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1A-7A.ABgenome.depth_4955sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/ab/chrAsub.ABgenome.depthVSsd.addMultipleColumn.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrAsub/chrAsub_AB_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1A-7A.ABgenome.depth_4955sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/ab/chrAsub.ABgenome.depthVSsd.addMultipleColumn.bin100_0.85.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrBsub/chrBsub_AB_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1B-7B.ABgenome.depth_4908sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/ab/chrBsub.ABgenome.depthVSsd.addMultipleColumn.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrBsub/chrBsub_AB_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1B-7B.ABgenome.depth_4908sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/ab/chrBsub.ABgenome.depthVSsd.addMultipleColumn.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/ab/chrBsub/chrBsub_AB_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/ab/001_chr1B-7B.ABgenome.depth_4908sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/ab/chrBsub.ABgenome.depthVSsd.addMultipleColumn.bin100_0.85.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/d/chrDgenome_bin100_0.65.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/d/002_chr1D-7D.Dgenome.depth_5049sites.txt", "6", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/d/chrDgenome.depthVSsd.addMultipleColumn.bin100_0.65.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/d/chrDgenome_bin100_0.75.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/d/002_chr1D-7D.Dgenome.depth_5049sites.txt", "7", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/d/chrDgenome.depthVSsd.addMultipleColumn.bin100_0.75.txt");
//        this.addGroupforsubsetDepthDB("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/003_getHighCumulativePos/d/chrDgenome_bin100_0.85.depthVSsd.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/004_subset5000sitefromDepthDB/d/002_chr1D-7D.Dgenome.depth_5049sites.txt", "8", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/008_cellMethod/005_subset.addMultipleColumn/d/chrDgenome.depthVSsd.addMultipleColumn.bin100_0.85.txt");
////
//        this.mergePosList("", "", "");
        //this.scriptReINFO();
        //this.bgzip_ABD();
        //new SplitScript().splitBwaScript("/Users/Aoyue/Documents/sh_bgzip_Vmap2_20191021.sh", "sh_bgzip_Vmap2_", 6 ,6 );
        //对过滤后的VCF文件进行质控,一、深度计算
        //new CountSites().calVcfAverageDepth("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/004_merged", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/006_depth_maf0.01");
        // 二、杂合度和缺失率计算
        //new CountSites().calSNPHetMissMaf("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/004_merged", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/007_mafHeterMiss_maf0.01");
        //new CountSites().calIndiHeter("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/005_all", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/007_mafHeterMiss_maf0.01");
        //this.extractMAF_D();
        //this.scriptFilterMafbyPopHexaDi();
//        this.scriptFilterMafbyPopHexaTetra();
        //this.extractMAF_D("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/004_maf001/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/002_calMAF/");
//        this.extractMAFbyPop_ABorD("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/005_maf001byPop", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/002_calMAF");
        //this.extractMAF("/Users/Aoyue/Documents/002_mergedbySub/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/009_calMAF_MAF0.01byPop");
//    new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/009_calMAF_MAF0.01bySub", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/010_bintable_bySub", "25", "0.5");
//        new CountSites().mergeTxt("/Users/Aoyue/Documents/test", "/Users/Aoyue/Documents/out/hexa.maf.txt"); //目的是将所有的六倍体MAF和在一起
//    new Bin().mkBarplotofMAF_single("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/012_calMAF_byPop/hexa.maf.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/013_bintable_byPop", "10", "0.5"); //hexa
//    new Bin().mkBarplotofMAF_single("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/012_calMAF_byPop/tetra.maf.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/013_bintable_byPop", "10", "0.5");
//    new Bin().mkBarplotofMAF_single("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/012_calMAF_byPop/di.maf.txt.gz", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/023_rebackDDtauschii/002_subsetVCFandMAF/013_bintable_byPop", "10", "0.5");
        //对新版的VMAPII进行计数
        //this.mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/012_countsiteFrommergedVCF/zzzzlog文件/log_021", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/012_countsiteFrommergedVCF/005_maf0.01SNPbyPop/001_filterIndelMaf0.01_20191016.txt");
//    new CountSites().mergeChr1and2txt_double("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/012_countsiteFrommergedVCF/005_maf0.01SNPbyPop/002_filterIndelMaf0.01_20191016.txt", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/012_countsiteFrommergedVCF/005_maf0.01SNPbyPop/003_filterIndelMaf0.01byRefChr_20191028.txt");
        /**
         * 对MAF>0.01byPo过滤的位点，进行最后的质控
         */
//    new CountSites().mergesubsetVCF("/Users/Aoyue/Documents/test/", "/Users/Aoyue/Documents/out/chr.hexa.vcf.gz");
//        this.filterNosegregation("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/009_databyPop/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/010_filterNosegregation/");
//this.filterNosegregation("/Users/Aoyue/Documents/test1/", "/Users/Aoyue/Documents/out1");
//    new CountSites().calSNPHetMissMaf("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/010_filterNosegregation", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/011_mafHeterMiss");
//        new CountSites().calxVcfAverageDepth("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/010_filterNosegregation", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/012_depth");
//        new CountSites().calIndiHeter("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/010_filterNosegregation", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/013_indiHeter");
//    new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/011_mafHeterMiss", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/014_bintable/001_maf", "10", "0.5");
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/011_mafHeterMiss", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/014_bintable/002_siteHeter", "50", "1");
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/011_mafHeterMiss", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/014_bintable/003_MissRate", "50", "1");
//            new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/012_depth", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/014_bintable/004_depth", "25", "25");
//        new Bin().mkBarplotofMAF("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/013_indiHeter", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/016_vcfQC/014_bintable/005_indiHeter", "50", "1");
//        this.scriptFilterMafbyPopHexaDi();
//        this.scriptFilterMafbyPopHexaTetra();
//        new SplitScript().splitBwaScript("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/014_filterVCF/013_fixVmap2/002_sh_filterMafbyPop_HexaTetraploid20191031.sh", "sh_filterMafbyPop_HexaTetraploid", 10, 3);
//        this.mergeTxt("/Users/Aoyue/Documents/log_023", "/Users/Aoyue/Documents/001.txt");
//        this.readVCF();
//        this.filterMissbyPopHexaDi("/Users/Aoyue/Documents/chr036.subgenome.maf0.01.SNP_bi.subset.vcf.gz", "/Users/Aoyue/Documents/chr036.subgenome.maf0.01.SNP_bi.subset_fileterMiss.vcf.gz");
//        this.filterMissbyPopHexaTetra("/Users/Aoyue/Documents/chr002.subgenome.maf0.01.SNP_bi.subset.vcf.gz", "/Users/Aoyue/Documents/chr001.ABgenome.filterMiss_subset.filterMiss.vcf.gz");
//    this.scriptFilterMiss();
//    new SplitScript().splitBwaScript("/Users/Aoyue/Desktop/sh_filterVCFbyMiss20191103.sh","sh_filterVCF",15,3);
//    this.mergeTxt("/Users/Aoyue/Documents/log_024", "/Users/Aoyue/Documents/lll.txt");
//        this.mergeVCFandFilter();


//        this.filterHeterbyPop();
//        this.runParallele_filterHeterbyPop();

    }

    public void runParallele_filterHeterbyPop(){ //本地运行常用
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/002_exonSNPVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/010_exonSNPVCF_filterHeter0.05";
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS, f.getName().split(".vcf")[0] + "_filterbyHeter0.05.vcf.gz").getAbsolutePath();
            this.filterHeterbyPop(infileS,outfileS);
            System.out.println(f.getName() + "\tis completed at " + outfileS);
        });
    }



    /**
     * 计算每个群体的缺失率，若都大于0.2，则过滤该位点
     *
     * @param infileS
     * @param outfileS
     */
    public void filterHeterbyPop(String infileS, String outfileS) {
        /**
         * 判断属于哪个基因组
         */
        String mapS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_chrList/ChrID.txt";
        String chr = new File(infileS).getName().substring(3,6);
        HashMap<String,String> hm = AoFile.getHashMapStringKey(mapS,1,5);
        String sub = hm.get(chr);
        String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/BreadWheat_S419.txt";
        String diFileS = null;

        if (sub.equals("AABB")){
            diFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/EmmerWheat_S187.txt";
        }
        else if (sub.equals("DD")){
            diFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/Ae.tauschii_S36.txt";
        }

//        String hexaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt";
//        String diFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/Ae.tauschii_S36.txt";

        String[] hexaArray = AoFile.getStringArraybyList_withoutHeader(hexaFileS,0);
        String[] diArray = AoFile.getStringArraybyList_withoutHeader(diFileS,0);
        List<Integer> indexHexa = new ArrayList<>();
        List<Integer> indexDi = new ArrayList<>();

        //预进行计算的数字：从过滤MAF0.01后的数据开始进行计算
//        System.out.println("Chr\tTotalSNP Num(MAF>0.01)\tBiallelic Num(MAF>0.01)\tTriallelic Num(MAF>0.01)\tIndel Num(MAF>0.01)\tInsertion Num(MAF>0.01)\tDeletion Num(MAF>0.01)");
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int ori = 0;
            int cntSNP = 0; //1.总共的SNP数量
            int cntBi = 0; //2.Bi-allelic SNP的数量
            int cntTri = 0; //3.Tr-allelic SNP的数量

            int cntBi_ABD = 0; //4.ABD六倍体的Bi-allelic SNP的数量
            int cntTri_ABD = 0; //5.ABD六倍体的Tri-allelic SNP的数量
            int cntBi_D = 0; //6.D二倍体的Bi-allelic SNP的数量
            int cntTri_D = 0; //7.D二倍体的Tri-allelic SNP的数量

            int cntIndel = 0; //8.总共Indel的数量
            int cntI = 0; //9.总共Insertion的数量
            int cntD = 0; //10.总共Deletion的数量

            int cntIndel_ABD = 0; //11.ABD六倍体的Indel的数量
            int cntI_ABD = 0; //12.ABD六倍体的Insertion的数量
            int cntD_ABD = 0; //13.ABD六倍体的Deletion的数量

            int cntIndel_D = 0; //14.D二倍体的Indel的数量
            int cntI_D = 0; //15.D二倍体的Insertion的数量
            int cntD_D = 0; //16.D二倍体的Deletion的数量

            while ((temp = br.readLine()) != null) {
                //***********************************************************//
                if (temp.startsWith("##")) {//将注释信息写入表格中
                    bw.write(temp);
                    bw.newLine();
                }

                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write(temp);
                    bw.newLine();
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(hexaArray, taxon);
                        int index2 = Arrays.binarySearch(diArray, taxon);

                        if (index1 > -1) {
                            indexHexa.add(i);
                        }
                        if (index2 > -1) {
                            indexDi.add(i);
                        }
                    }
                    Collections.sort(indexHexa);
                    Collections.sort(indexDi);
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    String alt = l.get(4);
                    List<String> lgeno = new ArrayList<>();
                    List<String> lHexaGeno = new ArrayList<>();
                    List<String> lDiGeno = new ArrayList<>();
                    String altList = l.get(4);
                    //先添加lgeno 再添加lHexaGeno lDiGeno
                    for (int i = 9; i < l.size(); i++) { //无论有无基因型，都加进去了
                        lgeno.add(l.get(i));
                    }
                    for (int i = 0; i < indexHexa.size(); i++) { //无论有无基因型，都加进去了
                        lHexaGeno.add(l.get(indexHexa.get(i)));
                    }
                    for (int i = 0; i < indexDi.size(); i++) { //无论有无基因型，都加进去了
                        lDiGeno.add(l.get(indexDi.get(i)));
                    }

                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    String[] hexaGenoArray = lHexaGeno.toArray(new String[lHexaGeno.size()]);
                    String[] diGenoArray = lDiGeno.toArray(new String[lDiGeno.size()]);

                    //如果该位点没有分离位点，则说明都是1/1 或都是 0/0，该群体在该位点的杂合度是1
//                    String missRate = this.getMissrate(genoArray);
//                    String missRateHexa = this.getMissrate(hexaGenoArray);
//                    String missRateDi = this.getMissrate(diGenoArray);

                    double heterRate = CalVCF.calSNPSitesHeter(genoArray);
                    double heterRateHexa = CalVCF.calSNPSitesHeter(hexaGenoArray);
                    double heterRateDi = CalVCF.calSNPSitesHeter(diGenoArray);

                    if (heterRateHexa < 0.05 && (heterRateDi < 0.05)) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(temp);
                        bw.write(sb.toString());
                        bw.newLine();

                        //保留任何一个群体中缺失率小于0.2的位点
                        if (!(alt.length() == 1)) { //2个alt的情况;若该位点含有D或I ，那么就属于Indel，如果没有D 或者I，那么就属于SNP
                            boolean ifD = false;
                            if (!altList.contains("D") && (!altList.contains("I"))) {
                                cntTri++;
                                cntSNP++;
                            }
                            if (altList.contains("D")) {
                                cntD++;
                                cntIndel++;
                                ifD = true;
                            }
                            if (altList.contains("I")) {
                                cntI++;
                                if (ifD == false) {
                                    cntIndel++;
                                }
                            }

                        } else if (alt.length() == 1) { //1个alt的情况;
                            if (!altList.equals("D") && (!altList.equals("I"))) {
                                cntBi++;
                                cntSNP++;
                            }
                            if (altList.equals("D")) {
                                cntD++;
                                cntIndel++;
                            }
                            if (altList.equals("I")) {
                                cntI++;
                                cntIndel++;
                            }
                        }
                    }
                } //
            }
            br.close();
            bw.flush();
            bw.close();
//            System.out.println(new File(infileS).getName().substring(3, 6) + "\t" + cntSNP + "\t" + cntBi + "\t" + cntTri + "\t" + cntIndel + "\t" + cntI + "\t" + cntD);
            System.out.println(infileS + " is completed at " + outfileS);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    
    
    
    

    public void scriptFilterMiss() {

        for (int i = 1; i < 43; i++) {
            int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
            Arrays.sort(db);
            if (Arrays.binarySearch(db, i) < 0) {
                String chr = PStringUtils.getNDigitNumber(3, i);//如果没有搜到，说明是不属于D的
                System.out.println("java -jar 024_filterMissbyPopHexaTetra.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/009_maf0.01SNPbyPopIncludeIndel/chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/010_miss0.2byPop/chr" + chr + "_miss0.2.vcf > log_024/log_filterMissbyPopHexaDi.chr" + chr + ".txt 2>&1");

            } else { //说明是属于D的
                String chr = PStringUtils.getNDigitNumber(3, i);
//                System.out.println("java -jar 021_filterMafbyPopHexaDi.jar /data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/006_bcftoolsMerge/chr" + chr + ".subgenome.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/008_maf0.01SNPbyPop/chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf > log_021/log_filterMafbyPopHexaDi.chr" + chr + ".txt");
                System.out.println("java -jar 024_filterMissbyPopHexaDi.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/009_maf0.01SNPbyPopIncludeIndel/chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/010_miss0.2byPop/chr" + chr + "_miss0.2.vcf > log_024/log_filterMissbyPopHexaDi.chr" + chr + ".txt 2>&1");

            }
        }
    }

    /**
     * 计算每个群体的缺失率，若都大于0.2，则过滤该位点
     *
     * @param infileS
     * @param outfileS
     */
    public void filterMissbyPopHexaTetra(String infileS, String outfileS) {

//        String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/BreadWheat_S419.txt";
//        String diFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/EmmerWheat_S187.txt";
        String hexaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt";
        String diFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/EmmerWheat_S187.txt";

        List<String> lhexa = new ArrayList<>(); //六倍体的taxa名集合
        List<String> ldi = new ArrayList<>(); //四倍体的taxa名集合
        List<Integer> indexHexa = new ArrayList<>();
        List<Integer> indexDi = new ArrayList<>();

        try {
            BufferedReader br = IOUtils.getTextReader(hexaFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                lhexa.add(temp);
            }
            br.close();
            br = IOUtils.getTextReader(diFileS);
            while ((temp = br.readLine()) != null) {
                ldi.add(temp);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        String[] hexaArray = lhexa.toArray(new String[lhexa.size()]);
        String[] diArray = ldi.toArray(new String[ldi.size()]);
        Arrays.sort(hexaArray);
        Arrays.sort(diArray);

        //预进行计算的数字：从过滤MAF0.01后的数据开始进行计算
        System.out.println("Chr\tTotalSNP Num(MAF>0.01)\tBiallelic Num(MAF>0.01)\tTriallelic Num(MAF>0.01)\tIndel Num(MAF>0.01)\tInsertion Num(MAF>0.01)\tDeletion Num(MAF>0.01)");
        try {
            BufferedReader br = null;
            BufferedWriter bw = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            if (outfileS.endsWith(".vcf")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".vcf.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            String temp = null;
            List<String> l = new ArrayList<>();

            int ori = 0;

            int cntSNP = 0; //1.总共的SNP数量
            int cntBi = 0; //2.Bi-allelic SNP的数量
            int cntTri = 0; //3.Tr-allelic SNP的数量

            int cntBi_ABD = 0; //4.ABD六倍体的Bi-allelic SNP的数量
            int cntTri_ABD = 0; //5.ABD六倍体的Tri-allelic SNP的数量
            int cntBi_D = 0; //6.D二倍体的Bi-allelic SNP的数量
            int cntTri_D = 0; //7.D二倍体的Tri-allelic SNP的数量

            int cntIndel = 0; //8.总共Indel的数量
            int cntI = 0; //9.总共Insertion的数量
            int cntD = 0; //10.总共Deletion的数量

            int cntIndel_ABD = 0; //11.ABD六倍体的Indel的数量
            int cntI_ABD = 0; //12.ABD六倍体的Insertion的数量
            int cntD_ABD = 0; //13.ABD六倍体的Deletion的数量

            int cntIndel_D = 0; //14.D二倍体的Indel的数量
            int cntI_D = 0; //15.D二倍体的Insertion的数量
            int cntD_D = 0; //16.D二倍体的Deletion的数量

            while ((temp = br.readLine()) != null) {
                //***********************************************************//
                if (temp.startsWith("##")) {//将注释信息写入表格中
                    bw.write(temp);
                    bw.newLine();
                }

                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write(temp);
                    bw.newLine();
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(hexaArray, taxon);
                        int index2 = Arrays.binarySearch(diArray, taxon);

                        if (index1 > -1) {
                            indexHexa.add(i);
                        }
                        if (index2 > -1) {
                            indexDi.add(i);
                        }
                    }
                    Collections.sort(indexHexa);
                    Collections.sort(indexDi);
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    String alt = l.get(4);
                    List<String> lgeno = new ArrayList<>();
                    List<String> lHexaGeno = new ArrayList<>();
                    List<String> lDiGeno = new ArrayList<>();
                    String altList = l.get(4);
                    //先添加lgeno 再添加lHexaGeno lDiGeno
                    for (int i = 9; i < l.size(); i++) { //无论有无基因型，都加进去了
                        lgeno.add(l.get(i));
                    }
                    for (int i = 0; i < indexHexa.size(); i++) { //无论有无基因型，都加进去了
                        lHexaGeno.add(l.get(indexHexa.get(i)));
                    }
                    for (int i = 0; i < indexDi.size(); i++) { //无论有无基因型，都加进去了
                        lDiGeno.add(l.get(indexDi.get(i)));
                    }

                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    String[] hexaGenoArray = lHexaGeno.toArray(new String[lHexaGeno.size()]);
                    String[] diGenoArray = lDiGeno.toArray(new String[lDiGeno.size()]);

                    String missRate = this.getMissrate(genoArray);
                    String missRateHexa = this.getMissrate(hexaGenoArray);
                    String missRateDi = this.getMissrate(diGenoArray);

                    double missrate = Double.parseDouble(missRate);
                    double missrateHexa = Double.parseDouble(missRateHexa);
                    double missrateDi = Double.parseDouble(missRateDi);

                    if (missrateHexa < 0.2 || (missrateDi < 0.2)) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(temp);
                        bw.write(sb.toString());
                        bw.newLine();

                        //保留任何一个群体中缺失率小于0.2的位点
                        if (!(alt.length() == 1)) { //2个alt的情况;若该位点含有D或I ，那么就属于Indel，如果没有D 或者I，那么就属于SNP
                            boolean ifD = false;
                            if (!altList.contains("D") && (!altList.contains("I"))) {
                                cntTri++;
                                cntSNP++;
                            }
                            if (altList.contains("D")) {
                                cntD++;
                                cntIndel++;
                                ifD = true;
                            }
                            if (altList.contains("I")) {
                                cntI++;
                                if (ifD == false) {
                                    cntIndel++;
                                }
                            }

                        } else if (alt.length() == 1) { //1个alt的情况;
                            if (!altList.equals("D") && (!altList.equals("I"))) {
                                cntBi++;
                                cntSNP++;
                            }
                            if (altList.equals("D")) {
                                cntD++;
                                cntIndel++;
                            }
                            if (altList.equals("I")) {
                                cntI++;
                                cntIndel++;
                            }
                        }
                    }
                } //
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(new File(infileS).getName().substring(3, 6) + "\t" + cntSNP + "\t" + cntBi + "\t" + cntTri + "\t" + cntIndel + "\t" + cntI + "\t" + cntD);
            System.out.println(infileS + " is completed at " + outfileS);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 计算每个群体的缺失率，若都大于0.2，则过滤该位点
     *
     * @param infileS
     * @param outfileS
     */
    public void filterMissbyPopHexaDi(String infileS, String outfileS) {
//        String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/BreadWheat_S419.txt";
//        String diFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/Ae.tauschii_S36.txt";

        String hexaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt";
        String diFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/Ae.tauschii_S36.txt";

        List<String> lhexa = new ArrayList<>(); //六倍体的taxa名集合
        List<String> ldi = new ArrayList<>(); //二倍体的taxa名集合
        List<Integer> indexHexa = new ArrayList<>();
        List<Integer> indexDi = new ArrayList<>();

        try {
            BufferedReader br = IOUtils.getTextReader(hexaFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                lhexa.add(temp);
            }
            br.close();
            br = IOUtils.getTextReader(diFileS);
            while ((temp = br.readLine()) != null) {
                ldi.add(temp);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        String[] hexaArray = lhexa.toArray(new String[lhexa.size()]);
        String[] diArray = ldi.toArray(new String[ldi.size()]);
        Arrays.sort(hexaArray);
        Arrays.sort(diArray);

        //预进行计算的数字：从过滤MAF0.01后的数据开始进行计算
        System.out.println("Chr\tTotalSNP Num(MAF>0.01)\tBiallelic Num(MAF>0.01)\tTriallelic Num(MAF>0.01)\tIndel Num(MAF>0.01)\tInsertion Num(MAF>0.01)\tDeletion Num(MAF>0.01)");
        try {
            BufferedReader br = null;
            BufferedWriter bw = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            if (outfileS.endsWith(".vcf")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".vcf.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            String temp = null;
            List<String> l = new ArrayList<>();

            int ori = 0;

            int cntSNP = 0; //1.总共的SNP数量
            int cntBi = 0; //2.Bi-allelic SNP的数量
            int cntTri = 0; //3.Tr-allelic SNP的数量

            int cntBi_ABD = 0; //4.ABD六倍体的Bi-allelic SNP的数量
            int cntTri_ABD = 0; //5.ABD六倍体的Tri-allelic SNP的数量
            int cntBi_D = 0; //6.D二倍体的Bi-allelic SNP的数量
            int cntTri_D = 0; //7.D二倍体的Tri-allelic SNP的数量

            int cntIndel = 0; //8.总共Indel的数量
            int cntI = 0; //9.总共Insertion的数量
            int cntD = 0; //10.总共Deletion的数量

            int cntIndel_ABD = 0; //11.ABD六倍体的Indel的数量
            int cntI_ABD = 0; //12.ABD六倍体的Insertion的数量
            int cntD_ABD = 0; //13.ABD六倍体的Deletion的数量

            int cntIndel_D = 0; //14.D二倍体的Indel的数量
            int cntI_D = 0; //15.D二倍体的Insertion的数量
            int cntD_D = 0; //16.D二倍体的Deletion的数量

            while ((temp = br.readLine()) != null) {
                //***********************************************************//
                if (temp.startsWith("##")) {//将注释信息写入表格中
                    bw.write(temp);
                    bw.newLine();
                }

                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write(temp);
                    bw.newLine();
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(hexaArray, taxon);
                        int index2 = Arrays.binarySearch(diArray, taxon);

                        if (index1 > -1) {
                            indexHexa.add(i);
                        }
                        if (index2 > -1) {
                            indexDi.add(i);
                        }
                    }
                    Collections.sort(indexHexa);
                    Collections.sort(indexDi);
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    String alt = l.get(4);
                    List<String> lgeno = new ArrayList<>();
                    List<String> lHexaGeno = new ArrayList<>();
                    List<String> lDiGeno = new ArrayList<>();
                    String altList = l.get(4);
                    //先添加lgeno 再添加lHexaGeno lDiGeno
                    for (int i = 9; i < l.size(); i++) { //无论有无基因型，都加进去了
                        lgeno.add(l.get(i));
                    }
                    for (int i = 0; i < indexHexa.size(); i++) { //无论有无基因型，都加进去了
                        lHexaGeno.add(l.get(indexHexa.get(i)));
                    }
                    for (int i = 0; i < indexDi.size(); i++) { //无论有无基因型，都加进去了
                        lDiGeno.add(l.get(indexDi.get(i)));
                    }

                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    String[] hexaGenoArray = lHexaGeno.toArray(new String[lHexaGeno.size()]);
                    String[] diGenoArray = lDiGeno.toArray(new String[lDiGeno.size()]);

                    String missRate = this.getMissrate(genoArray);
                    String missRateHexa = this.getMissrate(hexaGenoArray);
                    String missRateDi = this.getMissrate(diGenoArray);

                    double missrate = Double.parseDouble(missRate);
                    double missrateHexa = Double.parseDouble(missRateHexa);
                    double missrateDi = Double.parseDouble(missRateDi);

                    if (missrateHexa < 0.2 || (missrateDi < 0.2)) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(temp);
                        bw.write(sb.toString());
                        bw.newLine();

                        //保留任何一个群体中缺失率小于0.2的位点
                        if (!(alt.length() == 1)) { //2个alt的情况;若该位点含有D或I ，那么就属于Indel，如果没有D 或者I，那么就属于SNP
                            boolean ifD = false;
                            if (!altList.contains("D") && (!altList.contains("I"))) {
                                cntTri++;
                                cntSNP++;
                            }
                            if (altList.contains("D")) {
                                cntD++;
                                cntIndel++;
                                ifD = true;
                            }
                            if (altList.contains("I")) {
                                cntI++;
                                if (ifD == false) {
                                    cntIndel++;
                                }
                            }

                        } else if (alt.length() == 1) { //1个alt的情况;
                            if (!altList.equals("D") && (!altList.equals("I"))) {
                                cntBi++;
                                cntSNP++;
                            }
                            if (altList.equals("D")) {
                                cntD++;
                                cntIndel++;
                            }
                            if (altList.equals("I")) {
                                cntI++;
                                cntIndel++;
                            }
                        }
                    }
                } //
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(new File(infileS).getName().substring(3, 6) + "\t" + cntSNP + "\t" + cntBi + "\t" + cntTri + "\t" + cntIndel + "\t" + cntI + "\t" + cntD);
            System.out.println(infileS + " is completed at " + outfileS);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 获取VCF文件的基因型
     */
    public void readVCF() {
        try {
            String infileS = "/Users/Aoyue/Documents/chr036.vcf.gz";
            String outfileS = "/Users/Aoyue/Documents/missRate.txt";
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("MissingRate");
            bw.newLine();

            List<String> l = new ArrayList<>();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    String altList = l.get(4);
                    if (altList.contains("D") || altList.contains("I")) {
                        System.out.println("Alt allele is " + altList);
                    }
//                    List<String> lgeno = new ArrayList<>();
//                    for (int i = 9; i < l.size(); i++) { //无论有无基因型，都加进去了
//                        lgeno.add(l.get(i));
//                    }
//                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);

//                    bw.write(this.getMissrate(genoArray));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 获取变异位点的缺失率，该方法已经过验证，和TASSEL一致
     *
     * @param genoArray
     * @return
     */
    public String getMissrate(String[] genoArray) {
        int nz = 0; //有基因型的个体数
        for (int i = 0; i < genoArray.length; i++) {
            if (!genoArray[i].startsWith(".")) {
                nz++;
            }
        }
        double missRate = (double) (genoArray.length - nz) / genoArray.length;
        String missrate = String.format("%.4f", missRate);
        return missrate;
    }

    /**
     * 目的：过滤在群体中没有分离的位点,三等位位点不考虑
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void filterNosegregation(String infileDirS, String outfileDirS) {
        new File(outfileDirS).mkdirs();

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {

            try {
                String infileS = f.getAbsolutePath();
                BufferedReader br = null;
                BufferedWriter bw = null;
                String outfileS = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".removeNosegregationSites.vcf.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".removeNosegregationSites.vcf.gz")).getAbsolutePath();
                }
                bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                String te[] = null;
                List<String> l = new ArrayList<>();
                int cntSNP = 0; //totalSNP
                int cntkept = 0;
                System.out.println("File" + "\t" + "TotalSNP" + "\t" + "KeptSNP" + "\t" + "NosegregationSNP");
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        cntSNP++;
                        l = PStringUtils.fastSplit(temp);
                        String altList = l.get(4);
                        List<String> lgeno = new ArrayList<>();
                        if (altList.contains("D") || altList.contains("I")) {
                            continue;
                        }
                        if (l.get(4).length() > 1) {
                            continue;
                        }
                        for (int i = 9; i < l.size(); i++) {
                            lgeno.add(l.get(i));
                        }
                        String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                        boolean segregation = this.ifSegregation(genoArray, altList);
                        if (segregation == false) { //
                            continue;
                        }
                        cntkept++;
                        StringBuilder sb = new StringBuilder();
                        sb.append(temp);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }

                br.close();
                bw.flush();
                bw.close();
                System.out.println(new File(infileS).getName() + "\t" + cntSNP + "\t" + cntkept + "\t" + (cntSNP - cntkept));
                System.out.println(infileS + " is completed at " + outfileS);

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });
    }

    /**
     * 判断该群体在一个位点是否有分离，若有分离则返回true,若无分离则返回false,只针对有1个alt的情况
     *
     * @param genoArray
     * @param altList
     * @return
     */
    public boolean ifSegregation(String[] genoArray, String altList) {
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                continue;
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); // c是基因型第j个数值
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
        }
        //判断该群体在一个位点是否有分离，若有分离则返回true,若无分离则返回false
        boolean a = false;
        if ((acCnt[0] == 0 && acCnt[1] > 0) || (acCnt[0] > 0 && acCnt[1] == 0)) {
            a = false;
        } else {
            a = true;
        }
        return a;
    }

    /**
     * 目的：将log_21中所有txt文本的chr pos位点信息合并成一个文件，即按照亚群MAF0.01过滤的正式版本的sites。
     *
     * @param infileDirS
     * @param outfileS
     */
    public void mergelogTxt(String infileDirS, String outfileS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        //System.out.println("Chr\tSNP_Num");

        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            ///读表头，在第4行
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            br.readLine();
            br.readLine();
            br.readLine();
            bw.write(br.readLine()); //第四行是表头
            bw.newLine();

            //读正文部分
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String temp = br.readLine(); //read header
                //int chrint = Integer.parseInt(fs[i].getName().substring(3, 6));
                int cnt = 0;
                for (int j = 0; j < 7; j++) {  //每个log文件有7行，我们只要第5行的数据
                    temp = br.readLine();
                    if (j == 3) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(temp);
                        bw.write(sb.toString());
                        bw.newLine();

                    }
                }
                System.out.println(String.valueOf(fs[i].getName()) + "\t" + cnt);
            }
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 目的：从按照任一群体的MAF大于0.01中提取 MAF,MAF_ABD,MAF_AB,MAF_D，这里INFO信息是AAF_D； 表头按照：Maf
     * Maf_ABD Maf_ABorD
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void extractMAF(String infileDirS, String outfileDirS) {
        new File(outfileDirS).mkdirs();

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("vcf", "maf.txt.gz")).getAbsolutePath();
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("vcf.gz", "maf.txt.gz")).getAbsolutePath();
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            try {
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                bw.write("Maf\tMaf_ABD\tMaf_ABorD");
                bw.newLine();
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    //随机抽样 0.002 双等位位点
//                double r = Math.random();
//                    if (r > 0.002) {
//                        continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
//                    }
                    l = PStringUtils.fastSplit(temp.substring(0, 200));
                    //去除三等位位点
                    if (l.get(4).contains(",")) {
                        continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                    }
                    String INFO = l.get(7);
                    String maf = INFO.split(";")[6].split("=")[1];
                    String aafABD = INFO.split(";")[7].split("=")[1];
                    String aafABorD = INFO.split(";")[8].split("=")[1];
                    String mafABD = null;
                    String mafABorD = null;
                    mafABD = aafABD;
                    mafABorD = aafABorD;
                    if (Double.parseDouble(aafABD) > 0.5) {
                        mafABD = String.format("%.4f", 1 - Double.parseDouble(aafABD));
                    }
                    if (Double.parseDouble(aafABorD) > 0.5) {
                        mafABorD = String.format("%.4f", 1 - Double.parseDouble(aafABorD));
                    }
                    if (mafABD.equals("0.0000") || mafABD.equals("1.0000")) {
                        mafABD = "NA";
                    }
                    if (mafABorD.equals("0.0000") || mafABorD.equals("1.0000")) {
                        mafABorD = "NA";
                    }
                    StringBuilder sb = new StringBuilder();
                    sb.append(maf).append("\t").append(mafABD).append("\t").append(mafABorD);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();

                System.out.println(f.getName() + " is completed at" + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });
    }

    /**
     * 目的：从按照任一群体的MAF大于0.01中提取二倍体D的AAF值，这里INFO信息是AAF_D
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void extractMAFbyPop(String infileDirS, String outfileDirS) {
        new File(outfileDirS).mkdirs();

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("vcf", "maf.txt.gz")).getAbsolutePath();
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("vcf.gz", "maf.txt.gz")).getAbsolutePath();
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            try {
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                bw.write("Maf");
                bw.newLine();
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    //随机抽样 0.002 双等位位点 
//                double r = Math.random();
//                    if (r > 0.002) {
//                        continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
//                    }
                    l = PStringUtils.fastSplit(temp.substring(0, 200));
                    //去除三等位位点
                    if (l.get(4).contains(",")) {
                        continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                    }
                    String INFO = l.get(7);
                    String aafABorD = INFO.split(";")[8].split("=")[1];
                    double mafD = Double.MIN_VALUE;
                    if (aafABorD.equals("0.0000")) {
                        continue;
                    }
                    if (aafABorD.equals("1.0000")) {
                        continue;
                    }
                    mafD = Double.parseDouble(aafABorD);
                    if (Double.parseDouble(aafABorD) > 0.5) {
                        mafD = 1 - Double.parseDouble(aafABorD);
                    }
                    bw.write(String.format("%.4f", mafD));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 目的：从按照总群体的MAF大于0.01中提取D的MAF值，这里INFO信息是MAF_D
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void extractMAF_D(String infileDirS, String outfileDirS) {

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("vcf", "maf.txt.gz")).getAbsolutePath();
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
                outfileS = new File(outfileDirS, f.getName().replaceFirst("vcf.gz", "maf.txt.gz")).getAbsolutePath();
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            try {
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                bw.write("Maf");
                bw.newLine();
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    //随机抽样 0.002 双等位位点 
//                double r = Math.random();
//                    if (r > 0.002) {
//                        continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
//                    }

                    l = PStringUtils.fastSplit(temp.substring(0, 200));
                    //去除三等位位点
                    if (l.get(4).contains(",")) {
                        continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                    }
                    String INFO = l.get(7);
                    String mafD = INFO.split("MAF_D=")[1];
                    if (!mafD.equals("0.0000")) {
                        bw.write(mafD);
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    public void scriptFilterMafbyPopHexaTetra() {

        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
            Arrays.sort(db);
            if (Arrays.binarySearch(db, i) < 0) { //如果没有搜到，说明是不属于D的
//                System.out.println("java -jar 021_filterMafbyPopHexaTetra.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/004_rawMergedVCF_removeBadTaxa/chr" + chr + ".lineage.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/008_maf0.01SNPbyPop/chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf > log_021/log_filterMafbyPopHexaTetra.chr" + chr + ".txt");
                System.out.println("java -jar 023_filterMafbyPopHexaTetraKeepIndel.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/004_rawMergedVCF_removeBadTaxa/chr" + chr + ".lineage.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/009_maf0.01SNPbyPop/chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf > log_023/log_filterMafbyPopHexaTetra.chr" + chr + ".txt"); // 2>&1 &

            } else { //说明是属于D的
                //System.out.println("java -jar 021_filterMafbyPopHexaDi.jar /data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/006_bcftoolsMerge/chr" + chr + ".subgenome.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/008_maf0.01SNPbyPop/chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf > log_021/log_filterMafbyPopHexaDi.chr" + chr + ".txt");
            }
        }
    }

    //要解决的问题：1.修改header;2.修改info;3.只要有一个群体中MAF大于0.01，该位点保留
    public void filterMafbyPopHexaTetra(String infileS, String outfileS) {
        //String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/BreadWheat_S419.txt";
        //String tetraFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/EmmerWheat_S187.txt";

        String hexaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt";
        String tetraFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/EmmerWheat_S187.txt";
        List<String> lhexa = new ArrayList<>();
        List<String> ltetra = new ArrayList<>();
        List<Integer> indexHexa = new ArrayList<>();
        List<Integer> indexTera = new ArrayList<>();

        try {
            BufferedReader br = IOUtils.getTextReader(hexaFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                lhexa.add(temp);
            }
            br.close();
            br = IOUtils.getTextReader(tetraFileS);
            while ((temp = br.readLine()) != null) {
                ltetra.add(temp);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        String[] hexaArray = lhexa.toArray(new String[lhexa.size()]);
        String[] tetraArray = ltetra.toArray(new String[ltetra.size()]);
        Arrays.sort(hexaArray);
        Arrays.sort(tetraArray);
        //System.out.println("Chr\tTotalSNP Num\tBiallelic Num(Maf>0.01)\tTriallelic Num(Maf>0.01)\tTriallelic Num(Alt2F>0.01)\tTriallelic Ratio(Maf>0.01)\tTriallelic Ratio(Alt2F>0.01)");
        System.out.println("Chr\tTotalSNP Num\tBiallelic Num\tBiallelic Num(Maf_ABD>0.01)\tBiallelic Num(Maf_D>0.01)\tTriallelic Num\tDeletion Num\tInsertion Num\tIndel Num");
        try {
            BufferedReader br = null;
            BufferedWriter bw = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            if (outfileS.endsWith(".vcf")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".vcf.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            String temp = null;
            String te[] = null;
            List<String> l = new ArrayList<>();
            int cntSNP = 0; //totalSNP
            int cntBi = 0;
            int cntkept = 0;
            int cntTri = 0;
            int cntABD = 0;
            int cntAB = 0;
            int cntIndel = 0;
            int cntI = 0;
            int cntD = 0;

            while ((temp = br.readLine()) != null) {
                //***********************************************************//
                if (temp.startsWith("##")) {//过滤含有#的注释部分
                    if (temp.equals("##ALT=<ID=DEL,Description=\"Deletion\">")) {
                        temp = br.readLine();
                        bw.write("##INFO=<ID=AAF_ABD,Number=1,Type=Float,Description=\"Alternative allele frequency on hexaploid bread wheat\">");
                        bw.newLine();
                        bw.write("##INFO=<ID=AAF_AB,Number=1,Type=Float,Description=\"Alternative allele frequency on tetraploid emmer wheat\">");
                        bw.newLine();
                        //bw.write("##INFO=<ID=MAF_diploid,Number=1,Type=Float,Description=\"Minor allele frequency on diploid Aegilops tauschii\">");
                        //bw.newLine();
                        bw.write("##ALT=<ID=D,Description=\"Deletion\">");
                        bw.newLine();
                        bw.write("##ALT=<ID=I,Description=\"Insertion\">");
                        bw.newLine();
                        bw.write("##Species=Wheat");
                        bw.newLine();
                        bw.write("##ReferenceGenome=iwgsc_refseqv1.0");
                        bw.newLine();
                        bw.write("##VariantsMapVersion=\"vmap2\"");
                        bw.newLine();
                    }
                    if (temp.equals("##ALT=<ID=INS,Description=\"Insertion\">")) {
                        continue;
                    }
                    if (temp.startsWith("##bcftools")) {
                        continue;
                    }
                    if (temp.startsWith("##contig")) {
                        continue;
                    }
                    bw.write(temp);
                    bw.newLine();
                }

                //***********************************************************//
                //先替换taxa名字，在六倍体二倍体中，没有需要替换的taxa，故不进行任何taxa处理
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 0; i < 9; i++) {
                        bw.write(l.get(i) + "\t");

                    }
                    StringBuilder sb = new StringBuilder();
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        if (taxon.equals("PI428082")) {
                            taxon = "PI428082_1";
                        }
                        if (taxon.equals("PI466959_2")) {
                            taxon = "PI466959";
                        }

                        sb.append(taxon).append("\t");
                        int index1 = Arrays.binarySearch(hexaArray, taxon);
                        int index2 = Arrays.binarySearch(tetraArray, taxon);

                        if (index1 > -1) {
                            indexHexa.add(i);
                        }
                        if (index2 > -1) {
                            indexTera.add(i);
                        }
                    }
                    sb.deleteCharAt(sb.length() - 1);
                    bw.write(sb.toString());
                    bw.newLine();
                    Collections.sort(indexHexa);
                    Collections.sort(indexTera);

                }
                if (!temp.startsWith("#")) { //
                    l = PStringUtils.fastSplit(temp);
                    List<String> lgeno = new ArrayList<>();
                    List<String> lHexaGeno = new ArrayList<>();
                    List<String> lTetraGeno = new ArrayList<>();
                    String altList = l.get(4);
                    //过滤含有 D 或者 I 的位点
//                    if (altList.contains("D") || altList.contains("I")) {
//                        continue;
//                    }

                    //先添加lgeno 再添加lHexaGeno lTetraGeno
                    for (int i = 9; i < l.size(); i++) {
                        lgeno.add(l.get(i));
                    }
                    for (int i = 0; i < indexHexa.size(); i++) {
                        lHexaGeno.add(l.get(indexHexa.get(i)));
                    }
                    for (int i = 0; i < indexTera.size(); i++) {
                        lTetraGeno.add(l.get(indexTera.get(i)));
                    }

                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    String[] hexaGenoArray = lHexaGeno.toArray(new String[lHexaGeno.size()]);
                    String[] tetraGenoArray = lTetraGeno.toArray(new String[lTetraGeno.size()]);

//                    String INFO = this.getInfo(genoArray, altList);
//                    String hexaMAF = this.getSubgenomeMAF(hexaGenoArray, altList);
//                    String tetraMAF = this.getSubgenomeMAF(diGenoArray, altList);
                    String INFO = this.getInfo(genoArray, altList);
                    String hexaAAF = this.getSubgenomeAAF(hexaGenoArray, altList).split(",")[0];
                    String tetraAAF = this.getSubgenomeAAF(tetraGenoArray, altList).split(",")[0];

                    String hexaMAF = this.getSubgenomeAAF(hexaGenoArray, altList).split(",")[1];
                    String tetraMAF = this.getSubgenomeAAF(tetraGenoArray, altList).split(",")[1];

                    //开始进行,maf判断判断
                    double hexamaf = Double.parseDouble(hexaMAF);
                    double tetramaf = Double.parseDouble(tetraMAF);
                    if (hexamaf > 0.01 || (tetramaf > 0.01)) {
                        cntkept++;
                        boolean ifD = false;
                        if (altList.contains("D")) {
                            cntD++;
                            cntIndel++;
                            ifD = true;
                        }

                        if (altList.contains("I")) {
                            cntI++;
                            if (ifD == false) {
                                cntIndel++;
                            }
                        }

                        if (altList.length() > 1) {
                            cntTri++;
                        }
                        if (altList.length() == 1) {
                            cntBi++; //最终保留的二等位位点的个数
                            if (hexamaf > 0.01) {
                                cntABD++;
                            }
                            if (tetramaf > 0.01) {
                                cntAB++;
                            }
                        }
                        StringBuilder sb = new StringBuilder();
                        for (int i = 0; i < 7; i++) {
                            sb.append(l.get(i)).append("\t");
                        }
                        sb.append(INFO).append(";AAF_ABD=").append(hexaMAF).append(";AAF_AB=").append(tetraMAF).append("\tGT:AD:PL");
                        for (int i = 9; i < l.size(); i++) {
                            sb.append("\t").append(l.get(i));
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                } //
            }
            br.close();
            bw.flush();
            bw.close();
            cntSNP = cntkept - cntIndel;
            System.out.println(infileS.substring(infileS.indexOf("chr"), infileS.indexOf("chr") + 6) + "\t" + cntSNP + "\t" + cntBi + "\t" + cntABD + "\t" + cntD + "\t" + cntTri + "\t" + cntD + "\t" + cntI + "\t" + cntIndel);
            System.out.println(infileS + " is completed at " + outfileS);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void scriptFilterMafbyPopHexaDi() {

        for (int i = 1; i < 43; i++) {
            int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
            Arrays.sort(db);
            if (Arrays.binarySearch(db, i) < 0) { //如果没有搜到，说明是不属于D的
                continue;
            } else { //说明是属于D的
                String chr = PStringUtils.getNDigitNumber(3, i);
//                System.out.println("java -jar 021_filterMafbyPopHexaDi.jar /data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/006_bcftoolsMerge/chr" + chr + ".subgenome.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/008_maf0.01SNPbyPop/chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf > log_021/log_filterMafbyPopHexaDi.chr" + chr + ".txt");
//                System.out.println("java -jar 023_filterMafbyPopHexaDiKeepIndel.jar /data4/home/aoyue/vmap2/analysis/019_rebackDDtauschii/006_bcftoolsMerge/chr" + chr + ".subgenome.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/009_maf0.01SNPbyPop/chr" + chr + ".subgenome.maf0.01byPop.SNP.vcf > log_023/log_filterMafbyPopHexaDi.chr" + chr + ".txt 2>&1 &");

                System.out.println("java -jar 023_filterMafbyPopHexaDiKeepIndel.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/101_rawMergedVCF/chr" + chr + ".vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/102_MAF0.01/chr" + chr + "_maf0.01byPop.vcf > /data4/home/aoyue/vmap2/genotype/mergedVCF/log_102/log_filterMafbyPopHexaDi.chr" + chr + ".txt 2>&1 &");

            }
        }
    }

    //要解决的问题：1.修改header;2.修改info;3.只要有一个群体中MAF大于0.01，该位点保留
    //
    /**
     * 将DD Aegilops重新进行FastCall后，要解决的问题：1.REINFO,亚群的频率写成AAF； 2.保留任一群体大于0.01的位点
     *
     * @param infileS
     * @param outfileS
     */
    public void filterMafbyPopHexaDi(String infileS, String outfileS) {
        //String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/BreadWheat_S419.txt";
        //String diFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/Ae.tauschii_S36.txt";

        String hexaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt";
        String diFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/Ae.tauschii_S36.txt";

        List<String> lhexa = new ArrayList<>();
        List<String> ldi = new ArrayList<>();
        List<Integer> indexHexa = new ArrayList<>();
        List<Integer> indexDi = new ArrayList<>();

        try {
            BufferedReader br = IOUtils.getTextReader(hexaFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                lhexa.add(temp);
            }
            br.close();
            br = IOUtils.getTextReader(diFileS);
            while ((temp = br.readLine()) != null) {
                ldi.add(temp);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        String[] hexaArray = lhexa.toArray(new String[lhexa.size()]);
        String[] diArray = ldi.toArray(new String[ldi.size()]);
        Arrays.sort(hexaArray);
        Arrays.sort(diArray);
        //System.out.println("Chr\tTotalSNP Num\tBiallelic Num(Maf>0.01)\tTriallelic Num(Maf>0.01)\tTriallelic Num(Alt2F>0.01)\tTriallelic Ratio(Maf>0.01)\tTriallelic Ratio(Alt2F>0.01)");
        System.out.println("Chr\tTotalSNP Num\tBiallelic Num\tBiallelic Num(Maf_ABD>0.01)\tBiallelic Num(Maf_D>0.01)\tTriallelic Num\tDeletion Num\tInsertion Num\tIndel Num");
        try {
            BufferedReader br = null;
            BufferedWriter bw = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            if (outfileS.endsWith(".vcf")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".vcf.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            String temp = null;
            List<String> l = new ArrayList<>();
            int cntSNP = 0; //totalSNP
            int cntBi = 0;
            int cntkept = 0;
            int cntTri = 0;
            int cntABD = 0;

            int cnt = 0;
            int cntIndel = 0;
            int cntI = 0;
            int cntD = 0;
            while ((temp = br.readLine()) != null) {
                //***********************************************************//
                if (temp.startsWith("##")) {//过滤含有#的注释部分
                    if (temp.equals("##ALT=<ID=DEL,Description=\"Deletion\">")) {
                        temp = br.readLine();
                        bw.write("##INFO=<ID=AAF_ABD,Number=1,Type=Float,Description=\"Alternative allele frequency on hexaploid bread wheat\">");
                        bw.newLine();
                        bw.write("##INFO=<ID=AAF_D,Number=1,Type=Float,Description=\"Alternative allele frequency on diploid Aegilops tauschii\">");
                        bw.newLine();
                        bw.write("##ALT=<ID=D,Description=\"Deletion\">");
                        bw.newLine();
                        bw.write("##ALT=<ID=I,Description=\"Insertion\">");
                        bw.newLine();
                        bw.write("##Species=Wheat");
                        bw.newLine();
                        bw.write("##ReferenceGenome=iwgsc_refseqv1.0");
                        bw.newLine();
                        bw.write("##VariantsMapVersion=\"vmap2\"");
                        bw.newLine();
                    }
                    if (temp.equals("##ALT=<ID=INS,Description=\"Insertion\">")) {
                        continue;
                    }
                    if (temp.startsWith("##bcftools")) {
                        continue;
                    }
                    if (temp.startsWith("##contig")) {
                        continue;
                    }
                    bw.write(temp);
                    bw.newLine();
                }

                //***********************************************************//
                //先替换taxa名字，在六倍体二倍体中，没有需要替换的taxa，故不进行任何taxa处理
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write(temp);
                    bw.newLine();
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(hexaArray, taxon);
                        int index2 = Arrays.binarySearch(diArray, taxon);

                        if (index1 > -1) {
                            indexHexa.add(i);
                        }
                        if (index2 > -1) {
                            indexDi.add(i);
                        }
                    }
                    Collections.sort(indexHexa);
                    Collections.sort(indexDi);
                }
                if (!temp.startsWith("#")) { //
                    l = PStringUtils.fastSplit(temp);
                    List<String> lgeno = new ArrayList<>();
                    List<String> lHexaGeno = new ArrayList<>();
                    List<String> lDiGeno = new ArrayList<>();
                    String altList = l.get(4);
                    //过滤含有 D 或者 I 的位点
//                    if (altList.contains("D") || altList.contains("I")) {
//                        continue;
//                    }

                    //先添加lgeno 再添加lHexaGeno lDiGeno
                    for (int i = 9; i < l.size(); i++) {
                        lgeno.add(l.get(i));
                    }
                    for (int i = 0; i < indexHexa.size(); i++) {
                        lHexaGeno.add(l.get(indexHexa.get(i)));
                    }
                    for (int i = 0; i < indexDi.size(); i++) {
                        lDiGeno.add(l.get(indexDi.get(i)));
                    }

                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    String[] hexaGenoArray = lHexaGeno.toArray(new String[lHexaGeno.size()]);
                    String[] diGenoArray = lDiGeno.toArray(new String[lDiGeno.size()]);

                    String INFO = this.getInfo(genoArray, altList);
                    String hexaAAF = this.getSubgenomeAAF(hexaGenoArray, altList).split(",")[0];
                    String diAAF = this.getSubgenomeAAF(diGenoArray, altList).split(",")[0];

                    String hexaMAF = this.getSubgenomeAAF(hexaGenoArray, altList).split(",")[1];
                    String diMAF = this.getSubgenomeAAF(diGenoArray, altList).split(",")[1];

                    //开始进行,maf判断
                    double hexamaf = Double.parseDouble(hexaMAF);
                    double dimaf = Double.parseDouble(diMAF);
                    if (hexamaf > 0.01 || (dimaf > 0.01)) {

                        boolean ifD = false;
                        if (altList.contains("D")) {
                            cntD++;
                            cntIndel++;
                            ifD = true;
                        }

                        if (altList.contains("I")) {
                            cntI++;
                            if (ifD == false) {
                                cntIndel++;
                            }
                        }

                        cntkept++;

                        if (altList.length() > 1) {
                            cntTri++;
                        }
                        if (altList.length() == 1) {
                            cntBi++;
                            if (hexamaf > 0.01) {
                                cntABD++;
                            }
                            if (dimaf > 0.01) {
                                cntD++;
                            }
                        }
                        StringBuilder sb = new StringBuilder();
                        for (int i = 0; i < 7; i++) {
                            sb.append(l.get(i)).append("\t");
                        }
                        sb.append(INFO).append(";AAF_ABD=").append(hexaAAF).append(";AAF_D=").append(diAAF).append("\tGT:AD:PL");
                        for (int i = 9; i < l.size(); i++) {
                            sb.append("\t").append(l.get(i));
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                } //
            }
            br.close();
            bw.flush();
            bw.close();
            cntSNP = cntkept - cntIndel;
            System.out.println(infileS.substring(infileS.indexOf("chr"), infileS.indexOf("chr") + 6) + "\t" + cntSNP + "\t" + cntBi + "\t" + cntABD + "\t" + cntD + "\t" + cntTri + "\t" + cntD + "\t" + cntI + "\t" + cntIndel);
            System.out.println(infileS + " is completed at " + outfileS);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void bgzip_ABD() {
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/007_reINFO";
            String outfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/008_VMapII";
            System.out.println("bgzip -c -@ 20 " + new File(infileDirS, "chr" + chr + ".subgenome.maf0.01.SNP.vcf").getAbsolutePath() + " > " + new File(outfileDirS, "chr" + chr + ".vmap2.vcf.gz").getAbsolutePath() + " && tabix -p vcf " + new File(outfileDirS, "chr" + chr + ".vmap2.vcf.gz").getAbsolutePath() + " &");
        }
    }

    public void scriptReINFO() {
        List<Integer> l = new ArrayList<>();
        int j = 5;
        l.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            l.add(j);
        }

        int k = 6;
        l.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            l.add(k);
        }
        Collections.sort(l);

        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            int index = Collections.binarySearch(l, i);
            if (index < 0) {
                //System.out.println("java -jar 020_reINFOHexaTetraPloid.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/005_maf0.01SNP/chr" + chr + ".lineage.maf0.01.SNP.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/007_reINFO/chr" + chr + ".subgenome.maf0.01.SNP.vcf > log_020/reINFOHexaTetraPloid_chr" + chr + ".txt &");
                System.out.println("java -jar 020_reINFOHexaTetraPloid.jar /data4/home/aoyue/vmap2/analysis/018_subsetvcf/002_singleChr0.001/chr" + chr + ".lineage.maf0.01.SNP_bi.subset.vcf.gz /data4/home/aoyue/vmap2/analysis/018_subsetvcf/002_singleChr/chr" + chr + ".subgenome.maf0.01.SNP_bi.subset.vcf.gz > log_020/subset_reINFOHexaTetraPloid_chr" + chr + ".txt &");
            } else {
                //System.out.println("java -jar 020_reINFOHexaTetraPloid.jar /data4/home/aoyue/vmap2/genotype/mergedVCF/005_maf0.01SNP/chr" + chr + ".lineage.maf0.01.SNP.vcf /data4/home/aoyue/vmap2/genotype/mergedVCF/007_reINFO/chr" + chr + ".subgenome.maf0.01.SNP.vcf > log_020/reINFOHexaTetraPloid_chr" + chr + ".txt &");
                System.out.println("java -jar 020_reINFOHexaDiPloid.jar /data4/home/aoyue/vmap2/analysis/018_subsetvcf/002_singleChr0.001/chr" + chr + ".lineage.maf0.01.SNP_bi.subset.vcf.gz /data4/home/aoyue/vmap2/analysis/018_subsetvcf/002_singleChr/chr" + chr + ".subgenome.maf0.01.SNP_bi.subset.vcf.gz > log_020/subset_reINFOHexaTetraPloid_chr" + chr + ".txt &");
            }
        }
    }

    public void reINFOHexaDiPloid(String infileS, String outfileS) {
        ////// 先建立数据库，进行六倍体和四倍体的taxa整理并排序；
        //String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/BreadWheat_S419.txt";
        //String diFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/Ae.tauschii_S36.txt";

        String hexaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt";
        String diFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/Ae.tauschii_S36.txt";

        List<String> lhexa = new ArrayList<>();
        List<String> ldi = new ArrayList<>();
        List<Integer> indexHexa = new ArrayList<>();
        List<Integer> indexDi = new ArrayList<>();

        try {
            BufferedReader br = IOUtils.getTextReader(hexaFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                lhexa.add(temp);
            }
            br.close();
            br = IOUtils.getTextReader(diFileS);
            while ((temp = br.readLine()) != null) {
                ldi.add(temp);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        String[] hexaArray = lhexa.toArray(new String[lhexa.size()]);
        String[] diArray = ldi.toArray(new String[ldi.size()]);
        Arrays.sort(hexaArray);
        Arrays.sort(diArray);
        //System.out.println("Chr\tTotalSNP Num\tBiallelic Num(Maf>0.01)\tTriallelic Num(Maf>0.01)\tTriallelic Num(Alt2F>0.01)\tTriallelic Ratio(Maf>0.01)\tTriallelic Ratio(Alt2F>0.01)");
        try {
            BufferedReader br = null;
            BufferedWriter bw = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            if (outfileS.endsWith(".vcf")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".vcf.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            String temp = null;
            String te[] = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
//            int biallelicMafmoreNum = 0;
//            int cntCmaf12 = 0; //alt1 +alt2 大于0.005的个数
//            int cntCmaf2 = 0; //alt2 大于 0.005的个数
//            double ProportionofTriallelicMafmoreNum = Double.MIN_VALUE;
//            double ProportionofTriallelicAlt2moreNum = Double.MIN_VALUE;
            while ((temp = br.readLine()) != null) {
                //***********************************************************//
                if (temp.startsWith("##")) {//过滤含有#的注释部分
                    if (temp.equals("##ALT=<ID=DEL,Description=\"Deletion\">")) {
                        temp = br.readLine();
                        bw.write("##INFO=<ID=MAF_ABD,Number=1,Type=Float,Description=\"Minor allele frequency on hexaploid bread wheat\">");
                        bw.newLine();
                        bw.write("##INFO=<ID=MAF_D,Number=1,Type=Float,Description=\"Minor allele frequency on diploid Aegilops tauschii\">");
                        bw.newLine();
                        //bw.write("##INFO=<ID=MAF_diploid,Number=1,Type=Float,Description=\"Minor allele frequency on diploid Aegilops tauschii\">");
                        //bw.newLine();
                        bw.write("##ALT=<ID=D,Description=\"Deletion\">");
                        bw.newLine();
                        bw.write("##ALT=<ID=I,Description=\"Insertion\">");
                        bw.newLine();
                        bw.write("##Species=Wheat");
                        bw.newLine();
                        bw.write("##ReferenceGenome=iwgsc_refseqv1.0");
                        bw.newLine();
                    }
                    if (temp.equals("##ALT=<ID=INS,Description=\"Insertion\">")) {
                        continue;
                    }
                    if (temp.startsWith("##bcftools")) {
                        continue;
                    }
                    if (temp.startsWith("##contig")) {
                        continue;
                    }
                    bw.write(temp);
                    bw.newLine();
                }

                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write(temp);
                    bw.newLine();
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(hexaArray, taxon);
                        int index2 = Arrays.binarySearch(diArray, taxon);

                        if (index1 > -1) {
                            indexHexa.add(i);
                        }
                        if (index2 > -1) {
                            indexDi.add(i);
                        }
                    }
                    Collections.sort(indexHexa);
                    Collections.sort(indexDi);

                }
                if (!temp.startsWith("#")) { //
                    l = PStringUtils.fastSplit(temp);
                    List<String> lgeno = new ArrayList<>();
                    List<String> lHexaGeno = new ArrayList<>();
                    List<String> lDiGeno = new ArrayList<>();
                    String altList = l.get(4);
                    for (int i = 9; i < l.size(); i++) {
                        lgeno.add(l.get(i));
                        int index1 = Collections.binarySearch(indexHexa, i);
                        int index2 = Collections.binarySearch(indexDi, i);
                        if (index1 > -1) {
                            lHexaGeno.add(l.get(i));
                        }
                        if (index2 > -1) {
                            lDiGeno.add(l.get(i));
                        }
                    }
                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    String[] hexaGenoArray = lHexaGeno.toArray(new String[lHexaGeno.size()]);
                    String[] diGenoArray = lDiGeno.toArray(new String[lDiGeno.size()]);

                    String INFO = this.getInfo(genoArray, altList);
                    String hexaMAF = this.getSubgenomeMAF(hexaGenoArray, altList);
                    String diMAF = this.getSubgenomeMAF(diGenoArray, altList);
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < 7; i++) {
                        sb.append(l.get(i)).append("\t");
                    }
                    sb.append(INFO).append(";MAF_ABD=").append(hexaMAF).append(";MAF_D=").append(diMAF).append("\tGT:AD:PL");
                    for (int i = 9; i < l.size(); i++) {
                        sb.append("\t").append(l.get(i));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                } //
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void reINFOHexaTetraPloid(String infileS, String outfileS) {
        ////// 先建立数据库，进行六倍体和四倍体的taxa整理并排序；
//        String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/BreadWheat_S419.txt";
//        String tetraFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/removeBadTaxa/EmmerWheat_S187.txt";

        String hexaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/BreadWheat_S419.txt";
        String tetraFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/EmmerWheat_S187.txt";
        List<String> lhexa = new ArrayList<>();
        List<String> ltetra = new ArrayList<>();
        List<Integer> indexHexa = new ArrayList<>();
        List<Integer> indexTetra = new ArrayList<>();

        try {
            BufferedReader br = IOUtils.getTextReader(hexaFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                lhexa.add(temp);
            }
            br.close();
            br = IOUtils.getTextReader(tetraFileS);
            while ((temp = br.readLine()) != null) {
                ltetra.add(temp);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        String[] hexaArray = lhexa.toArray(new String[lhexa.size()]);
        String[] tetraArray = ltetra.toArray(new String[ltetra.size()]);
        Arrays.sort(hexaArray);
        Arrays.sort(tetraArray);

        //System.out.println("Chr\tTotalSNP Num\tBiallelic Num(Maf>0.01)\tTriallelic Num(Maf>0.01)\tTriallelic Num(Alt2F>0.01)\tTriallelic Ratio(Maf>0.01)\tTriallelic Ratio(Alt2F>0.01)");
        try {
            BufferedReader br = null;
            BufferedWriter bw = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            if (outfileS.endsWith(".vcf")) {
                bw = IOUtils.getTextWriter(outfileS);
            } else if (outfileS.endsWith(".vcf.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            String temp = null;
            String te[] = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
//            int biallelicMafmoreNum = 0;
//            int cntCmaf12 = 0; //alt1 +alt2 大于0.005的个数
//            int cntCmaf2 = 0; //alt2 大于 0.005的个数
//            double ProportionofTriallelicMafmoreNum = Double.MIN_VALUE;
//            double ProportionofTriallelicAlt2moreNum = Double.MIN_VALUE;
            while ((temp = br.readLine()) != null) {
                //***********************************************************//
                if (temp.startsWith("##")) {//过滤含有#的注释部分
                    if (temp.equals("##ALT=<ID=DEL,Description=\"Deletion\">")) {
                        temp = br.readLine();
                        bw.write("##INFO=<ID=MAF_ABD,Number=1,Type=Float,Description=\"Minor allele frequency on hexaploid bread wheat\">");
                        bw.newLine();
                        bw.write("##INFO=<ID=MAF_AB,Number=1,Type=Float,Description=\"Minor allele frequency on tetraploid emmer wheat\">");
                        bw.newLine();
                        //bw.write("##INFO=<ID=MAF_diploid,Number=1,Type=Float,Description=\"Minor allele frequency on diploid Aegilops tauschii\">");
                        //bw.newLine();
                        bw.write("##ALT=<ID=D,Description=\"Deletion\">");
                        bw.newLine();
                        bw.write("##ALT=<ID=I,Description=\"Insertion\">");
                        bw.newLine();
                        bw.write("##Species=Wheat");
                        bw.newLine();
                        bw.write("##ReferenceGenome=iwgsc_refseqv1.0");
                        bw.newLine();
                    }
                    if (temp.equals("##ALT=<ID=INS,Description=\"Insertion\">")) {
                        continue;
                    }
                    if (temp.startsWith("##bcftools")) {
                        continue;
                    }
                    if (temp.startsWith("##contig")) {
                        continue;
                    }
                    bw.write(temp);
                    bw.newLine();
                }

                //***********************************************************//
                //开始处理taxa的问题，先把所有taxa放入array中，记住在temp中的index
                if (temp.startsWith("#CHROM")) {
                    l = PStringUtils.fastSplit(temp);
                    bw.write(temp);
                    bw.newLine();
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(hexaArray, taxon);
                        int index2 = Arrays.binarySearch(tetraArray, taxon);

                        if (index1 > -1) {
                            indexHexa.add(i);
                        }
                        if (index2 > -1) {
                            indexTetra.add(i);
                        }
                    }
                    Collections.sort(indexHexa);
                    Collections.sort(indexTetra);

                }
                if (!temp.startsWith("#")) { //
                    l = PStringUtils.fastSplit(temp);
                    List<String> lgeno = new ArrayList<>();
                    List<String> lHexaGeno = new ArrayList<>();
                    List<String> lTetraGeno = new ArrayList<>();
                    String altList = l.get(4);
                    for (int i = 9; i < l.size(); i++) {
                        lgeno.add(l.get(i));
                        int index1 = Collections.binarySearch(indexHexa, i);
                        int index2 = Collections.binarySearch(indexTetra, i);
                        if (index1 > -1) {
                            lHexaGeno.add(l.get(i));
                        }
                        if (index2 > -1) {
                            lTetraGeno.add(l.get(i));
                        }
                    }
                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    String[] hexaGenoArray = lHexaGeno.toArray(new String[lHexaGeno.size()]);
                    String[] tetraGenoArray = lTetraGeno.toArray(new String[lTetraGeno.size()]);

                    String INFO = this.getInfo(genoArray, altList);
                    String hexaMAF = this.getSubgenomeMAF(hexaGenoArray, altList);
                    String tetraMAF = this.getSubgenomeMAF(tetraGenoArray, altList);
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < 7; i++) {
                        sb.append(l.get(i)).append("\t");
                    }
                    sb.append(INFO).append(";MAF_ABD=").append(hexaMAF).append(";MAF_AB=").append(tetraMAF).append("\tGT:AD:PL");
                    for (int i = 9; i < l.size(); i++) {
                        sb.append("\t").append(l.get(i));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                } //
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public String getSubgenomeAAF(String[] PopGenoArray, String altList) {
        int dp = 0;
        int nz = 0;
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] adCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt]; //GN到底代表什么？
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < PopGenoArray.length; i++) {
            if (PopGenoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(PopGenoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
//            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
//            for (int j = 0; j < temList.size(); j++) {
//                int c = Integer.parseInt(temList.get(j)); //c是第j个allele的深度值。注意AD的第一个是ref，第二个是次等位位点的深度，第三个是最小等位位点的深度
//                dp += c; //dp是总深度
//                adCnt[j] += c; //adCnt[j] 是第j个allele的深度值的总和，AD按照 ref alt1 alt2排序
//            }
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); // c是基因型第j个数值
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
//            int index1 = Integer.parseInt(temList.get(0)); //
//            int index2 = Integer.parseInt(temList.get(1));
//            gnCnt[index1][index2]++; //gnCnt[][]是二维数组，代表alt的个数的矩阵，比如有1个alt，则gnCnt[][]有gnCnt[0][0]  gnCnt[0][1] gnCnt[1][0] gnCnt[1][1]
//            if (index1 != index2) {
//                ht++;
//            }
        }
//        nz = PopGenoArray.length - nz;
        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }
        float aaf = (float) ((double) acCnt[1] / sum);

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%.4f", aaf)).append(",").append(String.format("%.4f", maf)); //.append(";MAF=")根据实际情况书写MAF_ABD MAF_AB MAF_D
        return sb.toString();
    }

    public String getSubgenomeMAF(String[] PopGenoArray, String altList) {
        int dp = 0;
        int nz = 0;
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] adCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt]; //GN到底代表什么？
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < PopGenoArray.length; i++) {
            if (PopGenoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(PopGenoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
//            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
//            for (int j = 0; j < temList.size(); j++) {
//                int c = Integer.parseInt(temList.get(j)); //c是第j个allele的深度值。注意AD的第一个是ref，第二个是次等位位点的深度，第三个是最小等位位点的深度
//                dp += c; //dp是总深度
//                adCnt[j] += c; //adCnt[j] 是第j个allele的深度值的总和，AD按照 ref alt1 alt2排序
//            }
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); // c是基因型第j个数值
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
//            int index1 = Integer.parseInt(temList.get(0)); //
//            int index2 = Integer.parseInt(temList.get(1));
//            gnCnt[index1][index2]++; //gnCnt[][]是二维数组，代表alt的个数的矩阵，比如有1个alt，则gnCnt[][]有gnCnt[0][0]  gnCnt[0][1] gnCnt[1][0] gnCnt[1][1]
//            if (index1 != index2) {
//                ht++;
//            }
        }
//        nz = PopGenoArray.length - nz;
        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%.4f", maf)); //.append(";MAF=")根据实际情况书写MAF_ABD MAF_AB MAF_D
        return sb.toString();
    }

    public String getInfo(String[] genoArray, String altList) {
        int dp = 0; //总深度
        int nz = 0; //有基因型的个体数
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] adCnt = new int[1 + nAlt]; //所有包括ref和alt的深度统计
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的基因型统计
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt]; //GN到底代表什么？
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型AD还有PL的集合

            //先计算深度
            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); //c是第j个allele的深度值。注意AD的第一个是ref，第二个是次等位位点的深度，第三个是最小等位位点的深度
                dp += c; //dp是总深度
                adCnt[j] += c; //adCnt[j] 是第j个allele的深度值的总和，AD按照 ref alt1 alt2排序
            }

            //再计算基因型
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) { //0/0:13,0:0,4,25
                int c = Integer.parseInt(temList.get(j)); // c是基因型0 1 2 其中的一个
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
            int index1 = Integer.parseInt(temList.get(0)); //0/0基因型的
            int index2 = Integer.parseInt(temList.get(1));
            gnCnt[index1][index2]++; //gnCnt[][]是二维数组，代表alt的个数的矩阵，比如有1个alt，则gnCnt[][]有gnCnt[0][0]  gnCnt[0][1] gnCnt[1][0] gnCnt[1][1]
            if (index1 != index2) {
                ht++;
            }
        }
        nz = genoArray.length - nz;
        int sum = 0; //所有allele的总数
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }

        StringBuilder sb = new StringBuilder();

        sb.append("DP=").append(dp).append(";NZ=").append(nz).append(";AD=");
        for (int i = 0; i < adCnt.length; i++) {
            sb.append(adCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length() - 1); //删除最后一个字符","号
        sb.append(";AC=");
        for (int i = 1; i < acCnt.length; i++) {
            sb.append(acCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.append(";GN=");
        for (int i = 0; i < gnCnt.length; i++) { //二维数组的长度是第一维的长度，这里是2（只有1个alt） 或者3 (有2个alt)
            for (int j = i + 1; j < gnCnt.length; j++) {
                sb.append(gnCnt[i][j]).append(",");
            }
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.append(";HT=").append(ht).append(";MAF=").append(String.format("%.4f", maf));
        return sb.toString();
    }

    /**
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void extractHapPosfromposAllele(String infileDirS, String outfileDirS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst("PosAllele.txt", "HapPos.txt.gz")).getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().replaceFirst("PosAllele.txt.gz", "HapPos.txt.gz")).getAbsolutePath();
                }

                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int cnt = 0;
                //bw.write("Chr\tPos\n");
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("Chr")) {
                        continue;
                    }
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    sb = new StringBuilder(l.get(0));
                    sb.append("\t").append(l.get(1));
                    bw.write(sb.toString());
                    bw.newLine();
                    if (cnt % 10000 == 0) {
                        System.out.println("Output " + String.valueOf(cnt) + " SNPs");
                    }
                    cnt++;
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt) + " SNPs output from " + f.getAbsolutePath());

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });
    }

    /**
     * 将chr001.ABDgenome 和 chr001.ABgenome 合并起来，生成一个合并后的文件
     */
    public void mergePosList(String inFileS1, String inFileS2, String outfileS) {
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
            BufferedReader br = IOUtils.getTextReader(inFileS1);
            String temp = null;
            while ((temp = br.readLine()).startsWith("##")) {
            };
            taxaNum1 = temp.split("\t").length - 9;
            String[] tem = null;
            int cnt1 = 0;
            while ((temp = br.readLine()) != null) {
                cnt1++;
                temp = temp.substring(0, 150);
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

            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
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
                        fre[index] = fre1[j + 1] * weight1; //复制fre1[j+1]到 fre中去， 权重1等于 在ABD中群体的个数 除以 在ABD中群体的个数加上AB群体的个数之和
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
                            fre[index] = fre2[j + 1] * weight2; //如果没搜到，说明在ABD群体中没有发现，只在群体2中发现；
                        } else {
                            fre[index] = fre[index] + fre2[j + 1] * weight2; //说明在ABD和AB群体中都有，2者相加。
                        }
                    }

                    int[] indices = PArrayUtils.getIndicesByDescendingValue(fre);
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
    /**
     * singleStream下面有parallelStream方法
     *
     * @param infileS
     * @param outfileS
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

    // 对抽样的vcf文件进行每个位点深度和sd进行统计和Pvalue获取，制成一个表格，chr pos averageDepth SD PValue
    // 注意表格不要以#开头，否则被注释，看不到
    /**
     * parallelStream上面有singleStream方法
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void statVcfDepth_SD_PValue(String infileDirS, String outfileDirS) {
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

            // chrAsub.ABDgenome.depthVSsd.addMultipleColumn.bin100_0.55.txt 文件名字
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
     * 对抽样的vcf文件进行每个位点每个taxa的深度统计，制成一个表格，chr pos averageDepth SD taxaDepth .....
     * 注意表格不要以#开头，否则被注释，看不到
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
