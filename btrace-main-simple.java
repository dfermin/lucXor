/* BTrace Script Template */
import com.sun.btrace.annotations.*;
import static com.sun.btrace.BTraceUtils.*;

@BTrace
public class TracingScript {
	@OnMethod(clazz="lucxor.LucXor", method="/.*/")
    public static void onMethod() {
        println("Hello from method");
    }
}