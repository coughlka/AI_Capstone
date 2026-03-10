# Configure the Azure provider
terraform {
  required_providers {
    azurerm = {
      source  = "hashicorp/azurerm"
      version = "~> 3.114.0"
    }
  }
  backend "azurerm" {
    resource_group_name  = "PSU-AI894-State-RG"
    storage_account_name = "psuai894storage"
    container_name       = "terraform-state"
    key                  = "terraform.tfstate"
  }
  required_version = ">= 1.1.0"
}

provider "azurerm" {
  features {}
}
