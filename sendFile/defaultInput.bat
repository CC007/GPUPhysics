////////////////////////////////////////////////////////////////////////////////
// JScript Section

try
{
    var WshShell = WScript.CreateObject('WScript.Shell');
    var Title = WScript.Arguments.Item(0);
    var Message = WScript.Arguments.Item(1);
    WshShell.AppActivate(Title);
    WshShell.SendKeys(Message);
    WScript.Quit(0);
}
catch(e)
{
    WScript.Echo(e);
    WScript.Quit(1);
}
WScript.Quit(2);