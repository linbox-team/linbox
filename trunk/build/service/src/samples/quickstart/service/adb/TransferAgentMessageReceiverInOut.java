
        /**
        * TransferAgentMessageReceiverInOut.java
        *
        * This file was auto-generated from WSDL
        * by the Apache Axis2 version: 1.2 Apr 27, 2007 (04:35:37 IST)
        */
        package samples.quickstart.service.adb;

        /**
        *  TransferAgentMessageReceiverInOut message receiver
        */

        public class TransferAgentMessageReceiverInOut extends org.apache.axis2.receivers.AbstractInOutSyncMessageReceiver{


        public void invokeBusinessLogic(org.apache.axis2.context.MessageContext msgContext, org.apache.axis2.context.MessageContext newMsgContext)
        throws org.apache.axis2.AxisFault{

        try {

        // get the implementation class for the Web Service
        Object obj = getTheImplementationObject(msgContext);

        TransferAgentSkeletonInterface skel = (TransferAgentSkeletonInterface)obj;
        //Out Envelop
        org.apache.axiom.soap.SOAPEnvelope envelope = null;
        //Find the axisOperation that has been set by the Dispatch phase.
        org.apache.axis2.description.AxisOperation op = msgContext.getOperationContext().getAxisOperation();
        if (op == null) {
        throw new org.apache.axis2.AxisFault("Operation is not located, if this is doclit style the SOAP-ACTION should specified via the SOAP Action to use the RawXMLProvider");
        }

        java.lang.String methodName;
        if(op.getName() != null & (methodName = org.apache.axis2.util.JavaUtils.xmlNameToJava(op.getName().getLocalPart())) != null){

        

            if("determinant".equals(methodName)){
                
                samples.quickstart.service.adb.xsd.DeterminantResponse determinantResponse17 = null;
                        samples.quickstart.service.adb.xsd.Determinant wrappedParam =
                                                             (samples.quickstart.service.adb.xsd.Determinant)fromOM(
                                    msgContext.getEnvelope().getBody().getFirstElement(),
                                    samples.quickstart.service.adb.xsd.Determinant.class,
                                    getEnvelopeNamespaces(msgContext.getEnvelope()));
                                                
                                               determinantResponse17 =
                                                   
                                                   
                                                         skel.determinant(wrappedParam)
                                                    ;
                                            
                                        envelope = toEnvelope(getSOAPFactory(msgContext), determinantResponse17, false);
                                    } else 

            if("rank".equals(methodName)){
                
                samples.quickstart.service.adb.xsd.RankResponse rankResponse19 = null;
                        samples.quickstart.service.adb.xsd.Rank wrappedParam =
                                                             (samples.quickstart.service.adb.xsd.Rank)fromOM(
                                    msgContext.getEnvelope().getBody().getFirstElement(),
                                    samples.quickstart.service.adb.xsd.Rank.class,
                                    getEnvelopeNamespaces(msgContext.getEnvelope()));
                                                
                                               rankResponse19 =
                                                   
                                                   
                                                         skel.rank(wrappedParam)
                                                    ;
                                            
                                        envelope = toEnvelope(getSOAPFactory(msgContext), rankResponse19, false);
                                    } else 

            if("trace".equals(methodName)){
                
                samples.quickstart.service.adb.xsd.TraceResponse traceResponse21 = null;
                        samples.quickstart.service.adb.xsd.Trace wrappedParam =
                                                             (samples.quickstart.service.adb.xsd.Trace)fromOM(
                                    msgContext.getEnvelope().getBody().getFirstElement(),
                                    samples.quickstart.service.adb.xsd.Trace.class,
                                    getEnvelopeNamespaces(msgContext.getEnvelope()));
                                                
                                               traceResponse21 =
                                                   
                                                   
                                                         skel.trace(wrappedParam)
                                                    ;
                                            
                                        envelope = toEnvelope(getSOAPFactory(msgContext), traceResponse21, false);
                                    } else 

            if("valence".equals(methodName)){
                
                samples.quickstart.service.adb.xsd.ValenceResponse valenceResponse23 = null;
                        samples.quickstart.service.adb.xsd.Valence wrappedParam =
                                                             (samples.quickstart.service.adb.xsd.Valence)fromOM(
                                    msgContext.getEnvelope().getBody().getFirstElement(),
                                    samples.quickstart.service.adb.xsd.Valence.class,
                                    getEnvelopeNamespaces(msgContext.getEnvelope()));
                                                
                                               valenceResponse23 =
                                                   
                                                   
                                                         skel.valence(wrappedParam)
                                                    ;
                                            
                                        envelope = toEnvelope(getSOAPFactory(msgContext), valenceResponse23, false);
                                    
            } else {
              throw new RuntimeException("method not found");
            }
        

        newMsgContext.setEnvelope(envelope);
        }
        }
        catch (Exception e) {
        throw org.apache.axis2.AxisFault.makeFault(e);
        }
        }
        
        //
            private  org.apache.axiom.om.OMElement  toOM(samples.quickstart.service.adb.xsd.Rank param, boolean optimizeContent){
            
                     return param.getOMElement(samples.quickstart.service.adb.xsd.Rank.MY_QNAME,
                                  org.apache.axiom.om.OMAbstractFactory.getOMFactory());
                    

            }
        
            private  org.apache.axiom.om.OMElement  toOM(samples.quickstart.service.adb.xsd.RankResponse param, boolean optimizeContent){
            
                     return param.getOMElement(samples.quickstart.service.adb.xsd.RankResponse.MY_QNAME,
                                  org.apache.axiom.om.OMAbstractFactory.getOMFactory());
                    

            }
        
            private  org.apache.axiom.om.OMElement  toOM(samples.quickstart.service.adb.xsd.Trace param, boolean optimizeContent){
            
                     return param.getOMElement(samples.quickstart.service.adb.xsd.Trace.MY_QNAME,
                                  org.apache.axiom.om.OMAbstractFactory.getOMFactory());
                    

            }
        
            private  org.apache.axiom.om.OMElement  toOM(samples.quickstart.service.adb.xsd.TraceResponse param, boolean optimizeContent){
            
                     return param.getOMElement(samples.quickstart.service.adb.xsd.TraceResponse.MY_QNAME,
                                  org.apache.axiom.om.OMAbstractFactory.getOMFactory());
                    

            }
        
            private  org.apache.axiom.om.OMElement  toOM(samples.quickstart.service.adb.xsd.Valence param, boolean optimizeContent){
            
                     return param.getOMElement(samples.quickstart.service.adb.xsd.Valence.MY_QNAME,
                                  org.apache.axiom.om.OMAbstractFactory.getOMFactory());
                    

            }
        
            private  org.apache.axiom.om.OMElement  toOM(samples.quickstart.service.adb.xsd.ValenceResponse param, boolean optimizeContent){
            
                     return param.getOMElement(samples.quickstart.service.adb.xsd.ValenceResponse.MY_QNAME,
                                  org.apache.axiom.om.OMAbstractFactory.getOMFactory());
                    

            }
        
            private  org.apache.axiom.om.OMElement  toOM(samples.quickstart.service.adb.xsd.Determinant param, boolean optimizeContent){
            
                     return param.getOMElement(samples.quickstart.service.adb.xsd.Determinant.MY_QNAME,
                                  org.apache.axiom.om.OMAbstractFactory.getOMFactory());
                    

            }
        
            private  org.apache.axiom.om.OMElement  toOM(samples.quickstart.service.adb.xsd.DeterminantResponse param, boolean optimizeContent){
            
                     return param.getOMElement(samples.quickstart.service.adb.xsd.DeterminantResponse.MY_QNAME,
                                  org.apache.axiom.om.OMAbstractFactory.getOMFactory());
                    

            }
        
                    private  org.apache.axiom.soap.SOAPEnvelope toEnvelope(org.apache.axiom.soap.SOAPFactory factory, samples.quickstart.service.adb.xsd.RankResponse param, boolean optimizeContent){
                      org.apache.axiom.soap.SOAPEnvelope emptyEnvelope = factory.getDefaultEnvelope();
                       
                                emptyEnvelope.getBody().addChild(param.getOMElement(samples.quickstart.service.adb.xsd.RankResponse.MY_QNAME,factory));
                            

                     return emptyEnvelope;
                    }
                    
                    private  org.apache.axiom.soap.SOAPEnvelope toEnvelope(org.apache.axiom.soap.SOAPFactory factory, samples.quickstart.service.adb.xsd.TraceResponse param, boolean optimizeContent){
                      org.apache.axiom.soap.SOAPEnvelope emptyEnvelope = factory.getDefaultEnvelope();
                       
                                emptyEnvelope.getBody().addChild(param.getOMElement(samples.quickstart.service.adb.xsd.TraceResponse.MY_QNAME,factory));
                            

                     return emptyEnvelope;
                    }
                    
                    private  org.apache.axiom.soap.SOAPEnvelope toEnvelope(org.apache.axiom.soap.SOAPFactory factory, samples.quickstart.service.adb.xsd.ValenceResponse param, boolean optimizeContent){
                      org.apache.axiom.soap.SOAPEnvelope emptyEnvelope = factory.getDefaultEnvelope();
                       
                                emptyEnvelope.getBody().addChild(param.getOMElement(samples.quickstart.service.adb.xsd.ValenceResponse.MY_QNAME,factory));
                            

                     return emptyEnvelope;
                    }
                    
                    private  org.apache.axiom.soap.SOAPEnvelope toEnvelope(org.apache.axiom.soap.SOAPFactory factory, samples.quickstart.service.adb.xsd.DeterminantResponse param, boolean optimizeContent){
                      org.apache.axiom.soap.SOAPEnvelope emptyEnvelope = factory.getDefaultEnvelope();
                       
                                emptyEnvelope.getBody().addChild(param.getOMElement(samples.quickstart.service.adb.xsd.DeterminantResponse.MY_QNAME,factory));
                            

                     return emptyEnvelope;
                    }
                    


        /**
        *  get the default envelope
        */
        private org.apache.axiom.soap.SOAPEnvelope toEnvelope(org.apache.axiom.soap.SOAPFactory factory){
        return factory.getDefaultEnvelope();
        }


        private  java.lang.Object fromOM(
        org.apache.axiom.om.OMElement param,
        java.lang.Class type,
        java.util.Map extraNamespaces){

        try {
        
                if (samples.quickstart.service.adb.xsd.Rank.class.equals(type)){
                
                           return samples.quickstart.service.adb.xsd.Rank.Factory.parse(param.getXMLStreamReaderWithoutCaching());
                    

                }
           
                if (samples.quickstart.service.adb.xsd.RankResponse.class.equals(type)){
                
                           return samples.quickstart.service.adb.xsd.RankResponse.Factory.parse(param.getXMLStreamReaderWithoutCaching());
                    

                }
           
                if (samples.quickstart.service.adb.xsd.Trace.class.equals(type)){
                
                           return samples.quickstart.service.adb.xsd.Trace.Factory.parse(param.getXMLStreamReaderWithoutCaching());
                    

                }
           
                if (samples.quickstart.service.adb.xsd.TraceResponse.class.equals(type)){
                
                           return samples.quickstart.service.adb.xsd.TraceResponse.Factory.parse(param.getXMLStreamReaderWithoutCaching());
                    

                }
           
                if (samples.quickstart.service.adb.xsd.Valence.class.equals(type)){
                
                           return samples.quickstart.service.adb.xsd.Valence.Factory.parse(param.getXMLStreamReaderWithoutCaching());
                    

                }
           
                if (samples.quickstart.service.adb.xsd.ValenceResponse.class.equals(type)){
                
                           return samples.quickstart.service.adb.xsd.ValenceResponse.Factory.parse(param.getXMLStreamReaderWithoutCaching());
                    

                }
           
                if (samples.quickstart.service.adb.xsd.Determinant.class.equals(type)){
                
                           return samples.quickstart.service.adb.xsd.Determinant.Factory.parse(param.getXMLStreamReaderWithoutCaching());
                    

                }
           
                if (samples.quickstart.service.adb.xsd.DeterminantResponse.class.equals(type)){
                
                           return samples.quickstart.service.adb.xsd.DeterminantResponse.Factory.parse(param.getXMLStreamReaderWithoutCaching());
                    

                }
           
        } catch (Exception e) {
        throw new RuntimeException(e);
        }
           return null;
        }



    

        /**
        *  A utility method that copies the namepaces from the SOAPEnvelope
        */
        private java.util.Map getEnvelopeNamespaces(org.apache.axiom.soap.SOAPEnvelope env){
        java.util.Map returnMap = new java.util.HashMap();
        java.util.Iterator namespaceIterator = env.getAllDeclaredNamespaces();
        while (namespaceIterator.hasNext()) {
        org.apache.axiom.om.OMNamespace ns = (org.apache.axiom.om.OMNamespace) namespaceIterator.next();
        returnMap.put(ns.getPrefix(),ns.getNamespaceURI());
        }
        return returnMap;
        }

        private org.apache.axis2.AxisFault createAxisFault(java.lang.Exception e) {
        org.apache.axis2.AxisFault f;
        Throwable cause = e.getCause();
        if (cause != null) {
            f = new org.apache.axis2.AxisFault(e.getMessage(), cause);
        } else {
            f = new org.apache.axis2.AxisFault(e.getMessage());
        }

        return f;
    }

        }//end of class
    