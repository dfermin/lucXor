/* BTrace Script Template */
import com.sun.btrace.annotations.*;
import static com.sun.btrace.BTraceUtils.*;

@BTrace
public class TracingScript {
	@OnMethod(clazz="lucxor.LucXor", method="main",
              location=@Location(
                    value = Kind.CALL,
                    clazz = "extra.HelloWorld",
                    method = "/call.*/",
                    where = Where.BEFORE)
    )
    public static void onMethod(
                  @ProbeMethodName(fqn = true) String pmn, 
                  @TargetMethodOrField(fqn = true) String tpmn) {
        println("Hello from method " + pmn);
        println("Going to invoke method " + tpmn);
    }
}