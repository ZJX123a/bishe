package test;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

public class testBaseOptions {
	/**
	 * 对baseOptions类中的基本操作进行测试
	 */
	baseOptions bo=new baseOptions();
	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void testContainNotAGCT() {
		String testKmer="ACGTAGCTGGGATTTCGAGAGGA";
		//System.out.println(bo.containNotAGCT(testKmer));
		assertEquals(false,bo.containNotAGCT(testKmer));
		testKmer="ACGTCGNGATATACTATATATTATA";
		assertEquals(true,bo.containNotAGCT(testKmer));
	}

	@Test
	public void testRevcomp() {
		String testKmer="ACTGCGATGC";
		//System.out.println(bo.revcomp(testKmer));
		assertEquals("GCATCGCAGT",bo.revcomp(testKmer));
		testKmer="CAGCATCAGGTCTCCAGAGCTGCAGAAGACGACGGCCGACTTGGATCACACTCTTGTGAG";
		assertEquals("CTCACAAGAGTGTGATCCAAGTCGGCCGTCGTCTTCTGCAGCTCTGGAGACCTGATGCTG",bo.revcomp(testKmer));
	}

	@Test
	public void testIntToBase() {
		assertEquals('A',bo.intToBase(0));
		assertEquals('C',bo.intToBase(1));
		assertEquals('G',bo.intToBase(2));
		assertEquals('T',bo.intToBase(3));
		assertEquals('Z',bo.intToBase(9));
		assertEquals('Z',bo.intToBase(-1));
		assertEquals('Z',bo.intToBase(4));
	}

	@Test
	public void testBaseToInt() {
		assertEquals(0,bo.baseToInt('A'));
		assertEquals(1,bo.baseToInt('C'));
		assertEquals(2,bo.baseToInt('G'));
		assertEquals(3,bo.baseToInt('T'));
		assertEquals(-1,bo.baseToInt('Q'));
		assertEquals(0,bo.baseToInt('a'));
		assertEquals(-1,bo.baseToInt('2'));
	//	fail("Not yet implemented");
	}

	@Test
	public void testIntvalToKmer() {
		assertEquals("ATCGCCT",bo.intvalToKmer(3479, 7));
		assertEquals("T",bo.intvalToKmer(3, 1));
	//	System.out.println(bo.kmerToIntval(""));
		assertEquals("ACTGCTGATCGATATAAAACGCCCA",bo.intvalToKmer(134017528628820l, 25));
	//	fail("Not yet implemented");
	} 

	@Test
	public void testComputeEntropy() {
		//System.out.println(bo.computeEntropy("C"));
		assertEquals(true,Math.abs(bo.computeEntropy("ATCGCCT")-1.842371)<0.0000001);
		assertEquals(true,Math.abs(bo.computeEntropy("ACTGCTGATCGATATAAAACGCCCA")-1.9322382)<0.0000001);
		assertEquals(true,Math.abs(bo.computeEntropy("C")-0.0)<0.0000001);
	//比较浮点数是否相等，可以比较其差是否小于某个界值
		//Math.abs(a-b)<0.00000001
	//	fail("Not yet implemented");
	}

	@Test
	public void testKmerToIntval() {
		//assertEquals("ATCGCCT",bo.intvalToKmer(3479, 7));
		assertEquals(3479,bo.kmerToIntval("ATCGCCT"));
		//assertEquals("T",bo.intvalToKmer(3, 1));
		assertEquals(3,bo.kmerToIntval("T"));
	//	System.out.println(bo.kmerToIntval(""));
	//	assertEquals("ACTGCTGATCGATATAAAACGCCCA",bo.intvalToKmer(134017528628820l, 25));
		assertEquals(134017528628820l,bo.kmerToIntval("ACTGCTGATCGATATAAAACGCCCA"));
	//	fail("Not yet implemented");
	}

}
