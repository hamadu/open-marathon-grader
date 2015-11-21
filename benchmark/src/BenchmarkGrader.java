import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;

public class BenchmarkGrader {
    public static void main(String[] args) {
        try {
            String exec = args[0];
            Long seed = Long.valueOf(args[1]);

            long limit = System.currentTimeMillis()+10000;
            Connector con = new Connector(exec);
            long ret = con.passIO();
            if (System.currentTimeMillis() > limit) {
                throw new RuntimeException("time up!");
            }
            System.out.println(ret);
        } catch (Exception e) {
            // something wrong. the score will be zero.
            System.out.println(0);
        }
    }

    public static class Connector {
        Runtime rt;
        Process proc;
        OutputStream os;
        BufferedReader br;

        public Connector(String exec) throws IOException {
            rt = Runtime.getRuntime();
            proc = rt.exec(exec);
            os = proc.getOutputStream();
            br = new BufferedReader(new InputStreamReader(proc.getInputStream()));
        }

        public long passIO() throws IOException {
            return Long.parseLong(br.readLine());
        }
    }
}
