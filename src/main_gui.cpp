
#include <QApplication>


#include "MainWindow.h"

int main(int argc, char** argv)
{
    
    QApplication app(argc, argv);
    app.setApplicationName("FiberCompression");
    MainWindow wnd;
    wnd.show();
    
    return app.exec();
    
            
}