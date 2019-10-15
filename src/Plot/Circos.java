/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Plot;

/**
 *
 * @author Aoyue
 */
public class Circos {

    public Circos() {
        this.band();

    }

    public void band() {
        for (int i = 1; i < 8; i++) {
            System.out.println("band\t" + i +  "A\t" + i + "AS\t" + i+ "AS\t0\t210200000\tgpos25");
            System.out.println("band\t" + i +  "B\t" + i + "BS\t" + i+ "BS\t0\t210200000\tgpos25");
            System.out.println("band\t" + i +  "D\t" + i + "DS\t" + i+ "DS\t0\t210200000\tgpos25");

        }

    }
}
