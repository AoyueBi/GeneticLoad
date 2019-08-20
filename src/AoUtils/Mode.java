/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class Mode {

    public Mode() {

    }

    public void mode() {
        String infileDirS = "";
        String outfileDirS = "";

        String infileS = "";
        String outfileS = "";

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                if (temp.startsWith("#")) {
                    continue;
                }
                cnt++;
            }
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        //////////////////
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        ////////////////////

    }

}
