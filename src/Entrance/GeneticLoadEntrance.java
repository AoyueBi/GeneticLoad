/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Entrance;

import WheatGeneticLoad.WheatABDvcfProcessor;
import AoUtils.SplitScript;
import WheatGeneticLoad.WheatABandDcleandataProcessor;
import WheatGeneticLoad.WheatBamDatabase;

/**
 *
 * @author Aoyue
 */
public class GeneticLoadEntrance {
    public GeneticLoadEntrance(){
        this.firstProcess();
        
        
    }
    public void firstProcess(){
        //new MapMake();
        //new Wheat120cleandataProcessor();
        //new Wheat120bamProcessor();
        //new Wheat200cleanDataProcessor();
        //new WheatBamDatabase();
        new WheatABandDcleandataProcessor();
        //new WheatABDvcfProcessor();
        
    }
    
    
    public static void main (String[] args){
        System.out.println("傲月-个人库 Here is the entrance of GeneticLoad!");
        //System.out.println("I made some revise on itellij ");
        new GeneticLoadEntrance();
        //new SplitScript().splitBwaScript("/Users/Aoyue/Documents/sh_md5_WheatVMapII_ABgenome_fixmatePosBam.sh", "md5_WheatVMapII_ABgenome_fixmateBam_", 20, 32);
        
     
        //new GeneticLoadEntrance().test("s", "s");
        
    }  
    
    public void test(String a, String b){
        //a = "aaa";
        //b = "bbb";
        System.out.println("a is " + a);
        System.out.println("b is " + b);
    }
}
