package svdpgm;
import java.io.IOException;
import svdpgm.SVDPGM;

public class Example {
	public static void main(String [] args) throws IOException
	{
		//Generate U, D and V transpose of example.pgm.
		SVDPGM.generateSVD("example.pgm");
		
		//Generate first 5 approximations of example.pgm using SVD
		for(int k = 1; k <= 5; k++){
			SVDPGM.generateApproximation("example.pgm", k);
		}
		
		//Generate pseudoinverse of exmple.pgm
		SVDPGM.generatePseudoinverse("example.pgm");
	}
}
