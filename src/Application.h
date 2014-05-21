
#ifndef APPLICATION_H
#define APPLICATION_H

#include <QCoreApplication>
#include <QString>
#include <QThread>

#include "../libs/tclap/CmdLine.h"

#include "Fibers.h"
#include "Utils.h"

class CmdLine;

class Application : public QCoreApplication
{
    Q_OBJECT
    
public:
    Application(int argc, char** argv);
    ~Application();

public slots:
    void showMessage(const QString&);
    void processFinished();
    void stopProcess();
    void close();
    
private:
    void connectSlots();
    bool initAndParseCommandLineOptions(int argc, char** argv);
    bool getAndValidateCommandLineArgs(Param& parameters);
    void startFromCommandLine(const Param& param);
    
    
    Fibers* m_pFibers;
    QThread* m_pThread;
    
    QFile         m_file;
    QElapsedTimer m_timer;
    
    Param         m_param;
    
    TCLAP::UnlabeledValueArg<std::string>* m_inputArg;
    TCLAP::UnlabeledValueArg<std::string>* m_outputArg;
    TCLAP::SwitchArg* m_forceArg;
    TCLAP::SwitchArg* m_compArg;
    TCLAP::SwitchArg* m_decompArg;
    TCLAP::ValueArg<float>* m_errorArg;
    TCLAP::ValueArg<int>* m_ttypeArg;
    TCLAP::ValueArg<int>* m_qtypeArg;
    TCLAP::ValueArg<int>* m_qnumArg;
    TCLAP::ValueArg<int>* m_etypeArg;
    TCLAP::ValueArg<std::string>* m_statsArg;
};

#endif // APPLICATION_H
