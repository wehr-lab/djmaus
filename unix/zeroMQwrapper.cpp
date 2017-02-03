/*
 * % Copyright (c) 2008 Shay Ohayon, California Institute of Technology.
 * % This file is a part of a free software. you can redistribute it and/or modify
 * % it under the terms of the GNU General Public License as published by
 * % the Free Software Foundation (see GPL.txt)
 */
//to compile with XCode 7.x on Mac you need to be running Matlab R2016a or higher mw070716

//note that mac/zeroMQwrapper.cpp and unix/zeroMQwrapper.cpp are identical
//but windows/zeroMQwrapper.cpp is slightly different for platform reasons

#include <stdio.h>
#include "mex.h"
#include "zmq.h"
#include "zmq_utils.h"
#include <string>
#include <cstring>
#include <queue>
#include <list>
#include <pthread.h>
#include <unistd.h>


typedef struct {
    const char *connect_url;
    unsigned long dwThread;
    bool thread_running;
    char *reply;
    
} ThreadData;

void* MyThreadFunction( void* lpParam );

#define MAX(x,y)(x>y)?(x):(y)
#define MIN(x,y)(x<y)?(x):(y)
#define ms100 100000
#define ms1000 1000000

std::queue<std::string> MyQueue;

typedef std::list<ThreadData*> lst;
lst MyHandleList;

bool initialized = false;
void* context;

void CloseContext(void)
{
    
    if (initialized) {
        
        for (lst::iterator i= MyHandleList.begin(); i != MyHandleList.end(); i++) {
            ThreadData* Tmp = *i;
            Tmp->thread_running = false;
            
        }
        usleep(ms100);
        zmq_ctx_destroy (context);
        initialized = false;
    }
}

void initialize_zmq()
{
    mexAtExit(CloseContext);
    context = zmq_ctx_new ();
    initialized = true;
}

void* MyThreadFunction( void* lpParam )
{
    ThreadData *pData = (ThreadData*) lpParam ;
    void *requester = zmq_socket (context, ZMQ_REQ);
    zmq_connect (requester, pData->connect_url); // "tcp://localhost:5555"
    pData->thread_running = true;
    
    while (pData->thread_running)
    {
        if (!MyQueue.empty())
        {
            
            std::string command = MyQueue.front();
            MyQueue.pop();
            
            zmq_msg_t request;
            zmq_msg_init_size (&request, command.length());
            memcpy (zmq_msg_data (&request), command.c_str() , command.length());
            zmq_msg_send (&request, requester, 0);
            
            zmq_msg_close (&request);
            
            zmq_msg_t reply;
            zmq_msg_init (&reply);
            zmq_msg_recv (&reply, requester, 0);
            
            //djmaus stuff:
            if   (strncmp(command.c_str(), "ChangeDirectory", 15) == 0)  {
                //if the message is "ChangeDirectory pathname"
                //this allows calling function to tell us to cd to pathname
                
                std::string mystr;
                mystr=command.c_str();
                char buffer[mystr.length()];
                std::size_t length = mystr.copy(buffer, mystr.length(), 16);
                buffer[length]='\0';
                if (chdir(buffer)== 0) {
                    printf("\nzmq wrapper: changed dir to %s", buffer);
                }
                else mexErrMsgIdAndTxt( "MATLAB:Could not change directory",
                        "Could not change directory.");
            }
            
            if   (strcmp(command.c_str(), "GetRecordingPath") == 0)  {
                //a way for calling function to learn the reply we get from GetRecordingPath
                
                // printf ("\nzmq wrapper: I will now attempt to write to RecordingPath.txt...");
                //    size_t mymsgsize=zmq_msg_size(&reply);
                //  char mymsg[mymsgsize]; bad idea
                //    memcpy (mymsg, zmq_msg_data (&reply) , mymsgsize);
                
//here's the new way to do it:
                size_t mymsgsize=zmq_msg_size(&reply);
                if (pData->reply == NULL) { //it should be empty
                    pData->reply = new char[mymsgsize + 1]; //+1 = room for termination
                    memcpy (pData->reply, zmq_msg_data (&reply) , mymsgsize);
                    pData->reply[mymsgsize] = '\0'; //terminate
                    
                    mexPrintf("\nzmq wrapper: GetRecordingPath returned %s \n", pData->reply);
                }
                else //if pData->reply is not NULL, something is wrong
                    mexErrMsgIdAndTxt("\nzmq wrapper: GetRecordingPath reply was not null??", "\nzmq wrapper: GetRecordingPath reply was not null??");
                
                //     FILE *fp;
                //    fp=fopen("RecordingPath.txt", "w");
                //   fprintf (fp, "%s", mymsg);
                //  fprintf (fp, "\n%d", mymsgsize);
                //  fclose(fp);
                //  printf ("\nzmq wrapper:wrote %s to RecordingPath.txt", mymsg);
                
                
            }
            
            zmq_msg_close (&reply);
            
        } else
        {
            // Be nice and sleep
            usleep(ms100); //was ms1000 mw070716
        }
    }
    
    zmq_close (requester);
    
    return 0;
    
}

ThreadData* create_client_thread(const char* connect_url, int n)
{
    
    ThreadData* pData = new ThreadData;
    pData->connect_url = new char[n];
    memcpy((char*)pData->connect_url,connect_url,n);
    pData->reply=NULL; //initialize to null
    pthread_t thread;
    
    int retValue = pthread_create(&thread,
            NULL,
            MyThreadFunction,
            (void*) pData);
    
    return pData;
    
}

void mexFunction( int nlhs, mxArray* plhs[],
        int nrhs, const mxArray* prhs[] )
{
    if (!initialized) initialize_zmq();
    
    if (nrhs < 2) {
        return;
    }
    
    
    char* Command = mxArrayToString(prhs[0]);
    
    if   (strcmp(Command, "StartConnectThread") == 0)  {
        int url_length = int(mxGetNumberOfElements(prhs[1])) + 1;
        char* url = mxArrayToString(prhs[1]);
        ThreadData* sd = create_client_thread(url, url_length);
        MyHandleList.push_back(sd);
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* Tmp= (double*)mxGetPr(plhs[0]);
        memcpy(Tmp, &sd, 8);
    }
    
    
    if (strcmp(Command, "Send") == 0)
    {
        double* Tmp= (double*) mxGetPr(prhs[1]);
        ThreadData* sd;
        memcpy(&sd, Tmp, 8);
        
        if (!sd->thread_running) {
            // Thread was killed and now we try to send things again?
            mexPrintf("Connection was closed. Cannot send.\n");
            return;
        }
        
        char *s= mxArrayToString(prhs[2]);
        std::string std_s(s);
        MyQueue.push(std_s);
        
        
    }
    
    
    if   (strcmp(Command, "GetReply") == 0)  { //get reply from ThreadData and put in plhs
        double* Tmp= (double*)mxGetPr(prhs[1]); //get handle
        ThreadData* sd;
        memcpy(&sd, Tmp, 8);
        if (sd->reply!=NULL) {//there is a reply waiting
            
            plhs[0] = mxCreateString(sd->reply); //put reply in plhs
            mexPrintf("\nzmq wrapper: mexFunction main putting %s in plhs\n", sd->reply);
            
            //now clear reply and deallocate
            delete[] sd->reply;
            sd->reply = NULL;
        }
        else //sd->reply is null so we cannot return a reply as requested
            mexPrintf("\nzmq wrapper GetReply: there is no reply available\n");
        
    }
    
    if (strcmp(Command, "CloseThread") == 0)
    {
        double* Tmp = (double*) mxGetPr(prhs[1]);
        ThreadData* sd;
        memcpy(&sd, Tmp, 8);
        sd->thread_running = false;
        
        // remove from active handle list
        MyHandleList.remove(sd);
    }
    
}
